# SEATREE modules: the research codes

This directory holds the compiled codes that the SEATREE modules
drive. They are separate programs with their own authors, licenses
and documentation; SEATREE's python drivers live in
`python3/seatree/modules/`. All of them are built from the top-level
`Makefile` (`make all`, `make optional`); see the main README.

Three of the codes are git submodules and are fetched with
`git submodule update --init` (or `git clone --recurse-submodules`).
The others are included in this repository.

## Mantle convection (`mc/`)

### hc (git submodule, pinned commit; upstream https://github.com/geodynamics/hc)

Global mantle circulation solver in the Hager and O'Connell (1981)
semi-analytical formulation: instantaneous flow and geoid for a given
density and radial viscosity structure, with free slip, no slip or
prescribed plate velocity surface boundary conditions. Written in C
by Thorsten Becker, based on Fortran routines by Brad Hager, Richard
O'Connell and Bernhard Steinberger. License: GPL. Citation and
references: `hc/README.md`. SEATREE builds it with `NO_GMT`; the GUI
uses `hc`, `hc_extract_sh_layer` and `sh_syn`.

### ConMan (git submodule, pinned commit; upstream https://github.com/geodynamics/conman)

Two-dimensional finite element code for incompressible, infinite
Prandtl number thermal convection, by Scott King, Arthur Raefsky and
Brad Hager; version 3.0 is distributed by the Computational
Infrastructure for Geodynamics. License: GPL v2 or later. Citation:
King, Raefsky and Hager (1990), *Phys. Earth Planet. Inter.*, 59,
195-207, and King, Raefsky and Hager (2020), ConMan v3.0.0,
doi:10.5281/zenodo.3633152. SEATREE generates input decks with
`ConMan/input_gen/make_conman_thermal_in.awk` and runs
`ConMan/conman`.

## Seismology (`seismo/`)

### syn2d

Two-dimensional Cartesian synthetic tomography for teaching: a
checkerboard (or an image) is used as input model, straight ray
synthetic travel times are generated with noise, and the model is
recovered by a damped least squares inversion. Fortran and C by Lapo
Boschi (2007) with modifications by Thorsten Becker (2008); see
`syn2d/readme.txt` and `doc/syn2d/`.

### larry

Surface wave phase velocity tomography in 2-D (pixel
parameterization, LSQR or Cholesky solvers, norm and roughness
damping). Fortran by Lapo Boschi; see `larry/README.TXT`. The phase
velocity measurements in `python3/data/larry/` are from Ekström et
al. (1997) and Boschi and Ekström (2002).

### larry3d

Three-dimensional body wave travel time tomography of the mantle
(voxel parameterization, LSQR with isotropic damping, after Boschi
and Dziewonski, 1999). Fortran by Lapo Boschi (2009). The travel
time data sets and reference models in `python3/data/larry3d/` are
described in that directory.

### nonlinloc

NonLinLoc, probabilistic earthquake location in 3-D and layered
models, by Anthony Lomax (copyright 1999-2007, GPL v2 or later).
Reference: Lomax et al. (2000), in *Advances in Seismic Event
Location*, Kluwer, 101-134; a copy and the user documentation are in
`nonlinloc/nll/doc/`. SEATREE provides a rudimentary interface for
locating events with the sample data in `python3/data/nll/`.

### EQdyna (git submodule, pinned commit; upstream https://github.com/EQDYNA/EQdyna)

Parallel finite element code for spontaneous dynamic rupture and
seismic wave propagation, by Benchun Duan, Dunyu Liu and co-workers
(Duan and Oglesby, 2006; Liu et al., 2020). License: MIT. The
submodule is pinned at the 5.3 era commit the SEATREE module was
written against; `../install/eqdyna_seatree_local_changes.patch`
holds the changes SEATREE needs on top of it until they are
upstream. Needs MPI (`mpif90`) and NetCDF Fortran; built with
`make eqdyna` (optional target). The module is an initial version.

### hypoinverse

Hypoinverse 2000 (version 1.1, 2007), the USGS earthquake location
program by Fred Klein (Klein, 2002, USGS Open-File Report 02-171).
Included for completeness and built with `make hypoinverse`
(optional target); there is no GUI module for it at present.

### tmp

A scratch copy of larry3d material from earlier development; not
used by the build or the GUI.
