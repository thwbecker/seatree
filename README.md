# SEATREE - Solid Earth Teaching and Research Environment

**Thorsten Becker** - thwbecker@post.harvard.edu - 11/2026

SEATREE is a modular, graphical front end to a set of solid Earth
research codes. It is meant for use in the classroom and for
interdisciplinary collaboration, and is used for one of the exercises
in *Tectonic Geodynamics* (Becker and Faccenna, Princeton University
Press, 2025).

The idea is to wrap existing, mostly Fortran and C, research codes in
transparent python drivers with a common graphical user interface,
so that a user can run them in a graphically supported "black box"
mode while still being able to see, and reuse, every command that was
executed. The research codes themselves are kept unchanged and
remain usable on their own; SEATREE adds the glue, the plotting, and
the GUI.

Project web site: https://twb.ig.utexas.edu/seatree/

## Contents

- [Modules and their authors](#modules-and-their-authors)
- [How SEATREE works](#how-seatree-works)
- [Requirements](#requirements)
- [Installation](#installation)
- [Running SEATREE](#running-seatree)
- [Build system](#build-system)
- [Troubleshooting](#troubleshooting)
- [Directory layout](#directory-layout)
- [History and credits](#history-and-credits)
- [How to cite](#how-to-cite)
- [License](#license)

## Modules and their authors

Each module wraps a separate code with its own authors, license and
citation. SEATREE only provides the interface. If you use results from
a module in a publication, please cite the underlying code as listed
in [How to cite](#how-to-cite) and in the module's own documentation.

| Module (GUI name) | Code and authors | Location | Build |
|---|---|---|---|
| HC-Flow Calculator: 3-D spherical mantle flow, geoid | **HC** by Thorsten Becker, after Hager and O'Connell (1981), with Fortran routines by Brad Hager, Richard O'Connell and Bernhard Steinberger | `modules/mc/hc`, git submodule of [geodynamics/hc](https://github.com/geodynamics/hc) | `make hc` |
| ConMan: 2-D thermal convection | **ConMan** by Scott King, Arthur Raefsky and Brad Hager (King et al., 1990; v3.0 2020) | `modules/mc/ConMan`, git submodule of [geodynamics/conman](https://github.com/geodynamics/conman) | `make conman` |
| Syn2D: 2-D Cartesian synthetic tomography (teaching) | **syn2d** by Lapo Boschi (2007), modified by Thorsten Becker (2008) | `modules/seismo/syn2d` | `make syn2d` |
| Larry: 2-D surface wave phase velocity tomography | **larry** by Lapo Boschi; phase velocity data from Ekström et al. (1997) and Boschi and Ekström (2002) | `modules/seismo/larry`, data in `python3/data/larry` | `make larry` |
| Larry 3-D: body wave mantle tomography | **larry3d** by Lapo Boschi (2009), after Boschi and Dziewonski (1999); travel time data sets as documented in `python3/data/larry3d` | `modules/seismo/larry3d` | `make larry3d` |
| NonLinLoc: earthquake location (rudimentary interface) | **NonLinLoc** by Anthony Lomax (Lomax et al., 2000) | `modules/seismo/nonlinloc/nll` | `make nll` |
| Earthquake Dynamics (EQDYNA): spontaneous rupture (initial version) | **EQdyna** by Benchun Duan, Dunyu Liu and co-workers (Duan and Oglesby, 2006; Liu et al., 2020) | `modules/seismo/EQdyna`, git submodule of [EQDYNA/EQdyna](https://github.com/EQDYNA/EQdyna) | `make eqdyna` (optional) |
| (no GUI module) | **Hypoinverse 2000** by Fred Klein, USGS (Klein, 2002) | `modules/seismo/hypoinverse` | `make hypoinverse` (optional) |

`make all` builds the first six, which the GUI needs. `make optional`
builds EQdyna, which needs an MPI Fortran compiler and the NetCDF
Fortran library, and hypoinverse, which is included for completeness.
See `modules/README.md` for licenses, versions and the origin of each
component.

## How SEATREE works

SEATREE is a python 3 / GTK 4 application. The parts are:

- **The main window** (`python3/seatree/gui/`): a start dialog to pick
  a module, a left panel with the module's controls, a plot area on
  the right, and a menu with Load Module, Save/Load Settings (as
  XML), Save Plot, and View Script.
- **Modules** (`python3/seatree/modules/<name>/`): one python package
  per research code. Each implements the small `Module` interface in
  `modules/module.py` (`setDefaults`, `getPanel`, `getPlotter`,
  `getSettings`, `loadSettings`, `cleanup`) and consists of a driver
  (for example `flowCalc.py`, `conman.py`), a GTK panel
  (`flowGUI.py`, `conmanGUI.py`), and, where needed, a thin wrapper
  that assembles the command lines of the underlying code. The module
  runs the compiled code as a subprocess, in a temporary directory
  under `/tmp`, on input files it writes from the GUI settings.
- **Plotters** (`python3/seatree/plotter/`): a GMT based plotter
  (used by HC, Larry, Larry 3-D, Syn2D, NonLinLoc) that writes
  PostScript with GMT and converts it to PNG for display, and a
  matplotlib based plotter (ConMan, EQdyna) that draws directly into
  the window.
- **The script log**: every shell command a module runs is recorded.
  "View Script" shows them, so a session can be turned into a shell
  script and repeated or modified outside the GUI. This is the main
  way the "black box" is kept transparent.
- **Configuration** (`python3/conf/`): `conf.xml` lists the modules
  and the GMT path, `<module>/<module>Conf.xml` holds the location of
  each module's binaries. Paths are relative to `python3/`. The
  installer writes these files; the committed defaults work for an
  unmoved checkout.
- **Data** (`python3/data/`): the input data each module offers in
  its menus (tomographic models, phase velocity measurements,
  viscosity and density profiles, plate boundaries, sample
  earthquake catalogs).

A new module is added by writing such a package, an installer class
(see `modules/example/`), and an entry in `python3/conf/install/`;
the Makefile gains a target for its compiled code.

## Requirements

- Linux (tested on Ubuntu 22.04 and 24.04) or macOS; Windows users
  can use WSL
- gcc, gfortran, make, git
- Python 3 with NumPy, SciPy, matplotlib (3.9 or newer), and the
  PyGObject bindings for GTK 4
- GMT 6 with the GSHHG coastlines, NetCDF, Ghostscript, and
  ImageMagick (`convert`) for the GMT based plots
- gawk, used to generate the ConMan input files
- optional, for EQdyna: MPICH (`mpif90`) and NetCDF Fortran
  (`libnetcdff-dev`)

HC is built without its own GMT bindings, so no GMT 4 development
files are needed. The older GMT 4 route is described below.

## Installation

### Quick start (Ubuntu)

```bash
# clone with the HC, ConMan and EQdyna submodules
git clone --recurse-submodules https://github.com/thwbecker/seatree.git
cd seatree

# packages (Ubuntu 22.04 or 24.04)
sudo apt-get install git make gcc gfortran gawk \
    python3 python3-numpy python3-scipy python3-matplotlib \
    python3-gi python3-gi-cairo gir1.2-gtk-4.0 \
    gmt gmt-dcw gmt-gshhg libgmt-dev libnetcdf-dev netcdf-bin \
    ghostscript imagemagick

# build HC, ConMan, syn2d, larry, larry3d and NonLinLoc
make all

# optional: EQdyna (needs mpich and libnetcdff-dev) and hypoinverse
make optional

# optional: HC regression tests and a short ConMan run
make test

# write the configuration files (and rebuild anything not built yet;
# answer the prompts, or accept the defaults with: yes '' | ...)
./install/configure.python3.gtk4

# check the installation
python3 verify.seatree.py
```

If the repository was cloned without `--recurse-submodules`, run
`git submodule update --init` first; `make` also does this when it
finds the submodule directories empty.

### macOS

With Homebrew: `brew install gcc gmt gtk4 pygobject3 ghostscript
imagemagick gawk netcdf` (and `mpich netcdf-fortran` for EQdyna),
`pip3 install matplotlib scipy`, then the same `make` and configure
steps as above. `./install.seatree.sh -e macos` runs the Homebrew
installation for you. macOS has seen less testing than Ubuntu.

### All in one script

`./install.seatree.sh -e ubuntu` (or `-e macos`) installs the
packages, runs `make all`, the python installer, `make eqdyna`, and
`verify.seatree.py`; `-m ubuntu` skips the package installation. The
script writes a log to `tmp/`. It also sets the environment variables
the optional GMT 4 route needs, and can be sourced for that purpose.

### GMT 4 (optional)

Older SEATREE versions required GMT 4. Setting `GMTVERSION=4` when
running `install.seatree.sh` builds GMT 4.5.18 and NetCDF from source
into the SEATREE directory (`install/install.gmt4.ubuntu22.sh`,
`install/install.netcdf.ubuntu22.sh`). This route is kept for
reference, is not needed for the current modules, and has seen less
testing than the GMT 6 default.

### Docker and packaged builds

The `Dockerfile` builds an Ubuntu 24.04 image with the steps above;
the GUI needs an X display forwarded into the container. For AppImage
and macOS application bundles see `build.sh` and
`distribute/README.md`. Older releases were also distributed as
VirtualBox Linux images, see
https://twb.ig.utexas.edu/seatree/install.html; those contain the old
python 2 version and are no longer maintained.

## Running SEATREE

```bash
./seatree              # or: python3 run.seatree.py
```

Both work from any directory. The start dialog lists the modules;
after choosing one, the left panel holds its controls and the right
side shows the plot. Load Module in the menu switches modules without
restarting. Settings can be saved to and loaded from XML files, plots
can be saved as PNG or PostScript, and View Script shows the commands
run so far.

A headless check of the HC and ConMan modules is in
`test/gui_smoke_test.py` (needs `xvfb-run`):

```bash
xvfb-run -a python3 test/gui_smoke_test.py hc
xvfb-run -a python3 test/gui_smoke_test.py conman
```

## Build system

The top-level `Makefile` builds every compiled component and is the
single place where compilers and flags are set:

- `make hc`: HC in `modules/mc/hc` with `NO_GMT` into
  `modules/mc/hc/bin`. The GUI only uses the ASCII output of the HC
  tools (`hc`, `hc_extract_sh_layer`, `sh_syn`) and converts it with
  the system GMT. `ARCH` and `HC_HOME` are cleared for this build so
  that it does not interfere with a separate HC installation.
- `make conman`: the CIG ConMan in `modules/mc/ConMan/src`; the
  binary is `modules/mc/ConMan/conman`, input decks are generated with
  `modules/mc/ConMan/input_gen/make_conman_thermal_in.awk`.
- `make syn2d`, `make larry`, `make larry3d`, `make nll`: the
  seismology modules with their own makefiles; the compilers and flags
  (`F77`, `F90`, `FFLAGS`, `F_EXT_SOURCE_FLAG`, `ARCH`) are passed in,
  so nothing needs to be edited in the module directories.
- `make eqdyna`: EQdyna from the submodule, with `MACHINE=ubuntu` or
  `macos` picked by `uname`; the binary is
  `modules/seismo/EQdyna/bin/eqdyna`.
- `make hypoinverse`: `modules/seismo/hypoinverse/bin/hyp`.
- `make test`: HC regression tests (both solvers) and a short ConMan
  run (`install/test.conman.sh`).
- `make clean`: removes the objects and binaries of all of the above.

The compilers (`CC`, `FC`, `F77`, `F90`) default to the `gcc` and
`gfortran` found in the path and can be set on the command line or in
the environment, for example `make FC=/usr/bin/gfortran all`; bare
names are resolved to full paths. `CFLAGS` and `FFLAGS` from the
environment are used if set; the flags the legacy sources need
(`-fcommon`, `-fallow-argument-mismatch`, `-ffixed-line-length-132`)
are added to them by the Makefile.

The python installer (`install/configure.python3.gtk4`, which runs
`python3/seatree/install/install.py`) calls the make targets above for
each module when asked to build it, and writes the configuration
files in `python3/conf/`. Build products are not tracked in git.

The EQdyna submodule is pinned to the upstream commit the previously
vendored copy was based on (844cefe, EQdyna 5.3). Changes that had
been made to that copy inside SEATREE (macOS NetCDF handling, a
surface output fix, per-panel rupture plots used by the GUI, a
single-rank test case) are kept in
`install/eqdyna_seatree_local_changes.patch` until they are in the
upstream repository; without them the EQdyna module compiles and runs
but does not produce the individual plot panels the GUI lists.
Current upstream EQdyna (5.25) has a restructured source tree, so
moving the pin forward also means adapting the module.

## Troubleshooting

- *"Permission denied" from make when running a compiler*, often in
  a submodule: GNU make 4.3 runs simple commands without a shell and
  stops at a directory that has the same name as the command in an
  earlier `PATH` entry (for instance with `.` in the path). The
  Makefile passes full compiler paths to avoid this; remove or rename
  such a directory to fix it for other builds, too.
- *"Rank mismatch in argument"* from gfortran in larry or syn2d: the
  legacy Fortran 77 needs `-fallow-argument-mismatch`; the Makefile
  adds it to whatever `FFLAGS` you have exported. If you call a
  module makefile directly, pass it yourself.
- *HC tests report missing output, `../bin` not found*: an exported
  `HC_HOME` or `ARCH` sent the HC build elsewhere; the Makefile clears
  both, make sure you build through it.
- *HC or ConMan not found by the GUI*: run `make all`, then
  `python3 verify.seatree.py`, which checks every binary and
  configuration path.
- *GMT errors or an empty plot area*: `conf.xml` points to the system
  GMT by default; if a configured path does not exist the loader falls
  back to the system path and says so on the terminal. GMT needs the
  GSHHG coastline data (`gmt-gshhg` on Ubuntu).

## Directory layout

- `python3/`: the python package (`seatree/`) with the GUI, module
  drivers, plotters, and the configuration (`conf/`) and data
  (`data/`) directories; see `python3/README.md`
- `modules/mc/`: mantle convection codes (HC and ConMan submodules)
- `modules/seismo/`: seismology codes (EQdyna is a submodule); see
  `modules/README.md`
- `install/`: helper scripts used by `install.seatree.sh`, the ConMan
  build test, and the pending EQdyna patch
- `test/`: headless GUI smoke test
- `distribute/`: AppImage and macOS bundle scripts
- `doc/`: the 2009 Eos paper, screenshots, syn2d notes
- `cider/`: handouts from the CIDER summer school exercises
- `legacy/`: the previous python 2 / GTK 2 version, kept for reference

## History and credits

The original framework (python 2, GTK 2) and the first module
drivers were written in 2008 and 2009 under the lead of Kevin Milner,
with help from others including Thorsten Becker, who initiated the
project, and the co-authors of Milner et al. (2009): Lapo Boschi,
J. Sain, Danijel Schorlemmer and H. Waterhouse. The 2024 to 2025
update to python 3 and GTK 4, including the matplotlib plotter, the
packaging scripts and the EQdyna module, was done by Dunyu Liu. In
2026, Thorsten Becker, assisted by Claude (Anthropic), reorganized
the build system (git submodules for HC, ConMan and EQdyna, the
top-level Makefile), fixed the remaining installation and GUI issues,
and revised this documentation.

The research codes behind the modules are the work of their authors,
listed in [Modules and their authors](#modules-and-their-authors) and
in `modules/README.md`.

Partial support for the original development was provided through
NSF-CAREER.

## How to cite

For SEATREE itself:

> Milner, K., Becker, T. W., Boschi, L., Sain, J., Schorlemmer, D.,
> and Waterhouse, H. (2009), The Solid Earth Research and Teaching
> Environment: a new software framework to share research tools in
> the classroom and across disciplines, *Eos Trans. AGU*, 90(12).
> https://twb.ig.utexas.edu/preprints/mbbssw08.pdf

For the codes behind the modules, please cite as recommended in their
own repositories; the key references are:

- HC: Hager, B. H., and O'Connell, R. J. (1981), A simple global
  model of plate dynamics and mantle convection, *J. Geophys. Res.*,
  86, 4843-4867; and the software citation in `modules/mc/hc/README.md`
  (Hager, O'Connell, Steinberger and Becker, HC, Computational
  Infrastructure for Geodynamics). The geoid and flow setup follows
  Becker, T. W., and Boschi, L. (2002), *Geochem. Geophys. Geosyst.*,
  3(1), 1003, doi:10.1029/2001GC000168.
- ConMan: King, S. D., Raefsky, A., and Hager, B. H. (1990), ConMan:
  vectorizing a finite element code for incompressible two-dimensional
  convection in the Earth's mantle, *Phys. Earth Planet. Inter.*, 59,
  195-207; software: King, S. D., Raefsky, A., and Hager, B. H.
  (2020), ConMan v3.0.0, doi:10.5281/zenodo.3633152.
- Larry and Larry 3-D: Boschi, L., and Dziewonski, A. M. (1999), High-
  and low-resolution images of the Earth's mantle: implications of
  different approaches to tomographic modeling, *J. Geophys. Res.*,
  104, 25567-25594. Surface wave data: Ekström, G., Tromp, J., and
  Larson, E. W. F. (1997), Measurements and global models of surface
  wave propagation, *J. Geophys. Res.*, 102, 8137-8157 (see
  `python3/data/larry/ekstrom_jgr97.pdf`), and Boschi, L., and
  Ekström, G. (2002), *J. Geophys. Res.*, 107, doi:10.1029/2000JB000059.
- NonLinLoc: Lomax, A., Virieux, J., Volant, P., and Berge-Thierry, C.
  (2000), Probabilistic earthquake location in 3D and layered models,
  in *Advances in Seismic Event Location*, Thurber, C. H., and
  Rabinowitz, N. (eds.), Kluwer, 101-134 (copy in
  `modules/seismo/nonlinloc/nll/doc`).
- EQdyna: Duan, B., and Oglesby, D. D. (2006), Heterogeneous fault
  stresses from previous earthquakes and the effect on dynamics of
  parallel strike-slip faults, *J. Geophys. Res.*, 111, B05309; see
  `modules/seismo/EQdyna/README.md` for the current software citation.
- Hypoinverse: Klein, F. W. (2002), User's guide to HYPOINVERSE-2000,
  a Fortran program to solve for earthquake locations and magnitudes,
  *USGS Open-File Report* 02-171.

## License

SEATREE, meaning the python framework, the GUI, the build system and
the syn2d, larry and larry3d drivers and data preparation, is
available under the GNU General Public License v2; see
[LICENSE.md](LICENSE.md). Copyright (c) Thorsten Becker and
contributors, 2008-2026.

The codes in `modules/` keep their own licenses and copyrights: HC
and ConMan are GPL (see their repositories), NonLinLoc is GPL
(copyright Anthony Lomax), EQdyna is MIT (see its repository),
Hypoinverse is USGS software (see its README files), and syn2d,
larry and larry3d were contributed by Lapo Boschi for distribution
with SEATREE. See `modules/README.md`.

## Support and contributing

If you want to use SEATREE in a classroom setting, we might be able
to offer some installation support and always welcome feedback. If
you would like to add a module, see [How SEATREE works](#how-seatree-works)
and `python3/README.md`, and let us know.
