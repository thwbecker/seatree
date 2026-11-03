# SEATREE - Solid Earth Teaching and Research Environment

**Thorsten Becker** - thwbecker@post.harvard.edu - 11/2026

SEATREE is a modular and user-friendly software to facilitate using
solid Earth research tools in the classroom and for interdisciplinary
scientific collaboration and used as one of the exercises in Tectonic
Geodynamics (Becker and Faccenna, Princeton University Press, 2025).

We use python wrappers and make use of modern software design
concepts, while remaining compatible with traditional scientific
coding. Our goals are to provide a fully contained, yet transparent
package that lets users operate in an easy, graphically supported
"black box" mode, while also allowing to look under the hood. In the
long run, we envision SEATREE to contribute to new ways of sharing
scientific research, and making (numerical) experiments truly
reproducible again.

## Modules

SEATREE is module based. The current version includes:

| Module | Code | Location | Build |
|---|---|---|---|
| 3-D spherical mantle flow | HC | `modules/mc/hc` (git submodule, [geodynamics/hc](https://github.com/geodynamics/hc)) | `make hc` |
| 2-D mantle (thermal) convection | ConMan | `modules/mc/ConMan` (git submodule, [geodynamics/conman](https://github.com/geodynamics/conman)) | `make conman` |
| 2-D synthetic tomography (teaching) | syn2d | `modules/seismo/syn2d` | python installer |
| Surface wave phase velocity tomography | larry | `modules/seismo/larry` | python installer |
| 3-D body wave mantle tomography | larry3d | `modules/seismo/larry3d` | python installer |
| Earthquake location (rudimentary) | NonLinLoc | `modules/seismo/nonlinloc` | python installer |
| Rupture dynamics (initial version) | EQdyna | `modules/seismo/EQdyna` | `install-eqdyna.sh` |

The main software design consists of transparent python wrappers that
drive the modules, including a GMT plotting tool, a VTK/Paraview 3-D
visualization interface, and a GTK 4 graphical user interface.

## Requirements

- Linux (tested on Ubuntu 22.04 and 24.04) or macOS; Windows users
  can use WSL
- C and Fortran compilers (gcc, gfortran), make, git
- Python 3 with NumPy, SciPy, matplotlib (3.9 or newer), and the
  PyGObject bindings for GTK 4
- GMT 6 (system installation) with the GSHHG coastlines, NetCDF,
  Ghostscript, and ImageMagick (`convert`) for plotting
- gawk, used to generate the ConMan input files

HC is built without its own GMT bindings, so no GMT 4 development
files are needed. The optional GMT 4 route is described below.

## Installation

### Quick start (Ubuntu)

```bash
# clone with the HC and ConMan submodules
git clone --recurse-submodules https://github.com/thwbecker/seatree.git
cd seatree

# packages (Ubuntu 22.04 or 24.04)
sudo apt-get install git make gcc gfortran gawk \
    python3 python3-numpy python3-scipy python3-matplotlib \
    python3-gi python3-gi-cairo gir1.2-gtk-4.0 \
    gmt gmt-dcw gmt-gshhg libgmt-dev libnetcdf-dev netcdf-bin \
    ghostscript imagemagick

# build HC and ConMan
make all

# optional: run the HC regression tests and a short ConMan run
make test

# build the seismology modules and write the configuration files
# (answer the prompts, or accept the defaults with: yes '' | ...)
./install/configure.python3.gtk4

# check the installation
python3 verify.seatree.py
```

If the repository was cloned without `--recurse-submodules`, run
`git submodule update --init` first; `make` will also do this when it
finds the submodule directories empty.

The combined script `./install.seatree.sh -e ubuntu` runs the package
installation and the steps above in one go, and `-m ubuntu` skips the
package installation. `./install.seatree.sh -e macos` does the
equivalent with Homebrew. These scripts also handle the optional
bundled NetCDF and GMT 4 builds, see the comments in the script.

### Running SEATREE

```bash
./seatree              # or: python3 run.seatree.py
```

Both work from any directory. The start dialog lets you select a
module; the module can be changed later from the menu.

### Build details

The top-level `Makefile` has these targets:

- `make hc`: builds HC in `modules/mc/hc` with `NO_GMT` into
  `modules/mc/hc/bin`. The GUI only uses the ASCII output of the HC
  tools (`hc`, `hc_extract_sh_layer`, `sh_syn`) and converts it with
  the system GMT.
- `make conman`: builds the CIG ConMan in `modules/mc/ConMan/src`;
  the binary is `modules/mc/ConMan/conman`. Input decks are generated
  by `modules/mc/ConMan/input_gen/make_conman_thermal_in.awk`.
- `make test`: HC regression tests (both solvers) and a short ConMan
  run (`install/test.conman.sh`).
- `make clean`: removes the objects and binaries of both codes.

The compilers default to the `gcc` (HC) and `gfortran` (ConMan) found
in the path, passed as full paths, and can be set on the command line,
for example `make FC=/usr/bin/gfortran all`. (GNU make 4.3 reports
"Permission denied" when running a compiler by its bare name if a
directory with that name exists in an earlier path entry; the full
path avoids this.)

The python installer (`install/configure.python3.gtk4`, which runs
`python3/seatree/install/install.py`) calls `make hc` and
`make conman` when asked to build these modules, builds the
seismology modules with their own makefiles, and writes the
configuration files in `python3/conf/`. Paths in these files are
relative to `python3/`, so the tree can be moved.

A headless GUI test for the HC and ConMan modules is in
`test/gui_smoke_test.py`; it needs `xvfb-run`:

```bash
xvfb-run -a python3 test/gui_smoke_test.py hc
xvfb-run -a python3 test/gui_smoke_test.py conman
```

### GMT 4 (optional)

Older SEATREE versions required GMT 4. Setting `GMTVERSION=4` when
running `install.seatree.sh` builds GMT 4.5.18 and NetCDF from source
into the SEATREE directory (`install/install.gmt4.ubuntu22.sh`,
`install/install.netcdf.ubuntu22.sh`). This route is kept for
reference and is not needed for the current modules; it has seen
less testing than the GMT 6 default.

### Docker

The `Dockerfile` builds an Ubuntu 24.04 image with the same steps as
above. The GUI needs an X display forwarded into the container.

### Other installation options

For AppImage and macOS application bundles, see `build.sh` and
`distribute/README.md`. For older download options, see:

http://www-udc.ig.utexas.edu/external/becker/seatree/wiki/download.html

## Directory layout

- `python3/`: the python package (`seatree/`) with the GUI, module
  wrappers, plotters, and the configuration files in `conf/`
- `modules/mc/`: mantle convection codes (HC, ConMan submodules)
- `modules/seismo/`: seismology codes
- `install/`: helper scripts used by `install.seatree.sh`
- `test/`: GUI smoke test
- `distribute/`: AppImage and macOS bundle scripts
- `doc/`, `cider/`: documentation, screenshots, and course handouts
- `legacy/`: the previous python 2 / GTK 2 version, kept for reference

## Resources

**Project Website:** http://www-udc.ig.utexas.edu/external/becker/seatree/

## Reference

Milner, K., Becker, T. W., Boschi, L., Sain, J., Schorlemmer, D. and H. Waterhouse: The Solid Earth Research and Teaching Environment: a new software framework to share research tools in the classroom and across disciplines. *Eos Trans. AGU*, 90, 12, 2009.

http://www-udc.ig.utexas.edu/external/becker/preprints/mbbssw08.pdf

If you use the HC or ConMan modules in publications, please also
cite these codes as described in their own repositories.

## License

SEATREE is freely available under the GNU General Public License v2. See [LICENSE.md](LICENSE.md) for the full license text.

## Support and Contributing

If you want to use SEATREE in a classroom setting, we might be able to
offer you some installation support and always welcome your
feedback. Also, if you like to add your own module to SEATREE, please
let us know; we might be able to provide some assistance. The original
code framework is mainly due to Kevin Milner, with a major 2025 update
by Dunyu Liu.

## Acknowledgments

Partial support was provided through NSF-CAREER.

## Third-Party Components and Copyright

SEATREE includes and may be distributed with the following third-party components:

- **HC**: mantle circulation code (Becker and co-workers, based on Hager and O'Connell), GPL, https://github.com/geodynamics/hc

- **ConMan**: 2-D finite element convection code by King, Raefsky and Hager, GPL, https://github.com/geodynamics/conman; cite as King, S. D., A. Raefsky, and B. H. Hager (2020), ConMan v3.0.0 [software], doi:10.5281/zenodo.3633152.

- **fstrack**: Set of routines to trace particles in mantle flow. Routines use precomputed velocities on GMT grd files. See: Becker et al. (GJI, 155, 696, 2003).

- **d-rex**: DREX Kaminski & Ribe texture routines by Eduard Kaminski and Jules Browaeys. See: Kaminski et al. (GJI, 157, 1, 2004) and Browaeys & Chevrot (GJI, 159, 667, 2004).

- **single_layer and multi_layer**: Shear wave splitting routines by Vera Schulte-Pelkum. See: Schulte-Pelkum & Blackman (GJI, 154, 166, 2003).

- **eispack**: EISPACK eigensystem routines

- **menke_splitting**: Bill Menke's cross-correlation routines

Copyright for these components remains with the original authors.

For the rest: Copyright (c) Thorsten Becker 2004-2026
