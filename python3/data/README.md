# Data used by the SEATREE modules

- `common/`: plate boundaries (NUVEL) for map overlays
- `hc/`: for the HC module: tomographic models as spherical harmonic
  expansions (`tomography/`, named after their publications), radial
  viscosity (`viscosity/`) and density scaling (`dscale/`) profiles,
  PREM, plate velocities (`pvelocity/`), and reference geoids
- `larry/`: Love and Rayleigh wave phase velocity measurements at
  several periods (Ekström et al., 1997; Boschi and Ekström, 2002,
  see the included paper) for the Larry module
- `larry3d/`: body wave travel time data sets for the Larry 3-D
  module, see `larry3d/README.md`
- `nll/`: sample velocity model and phase data for the NonLinLoc
  module
- `syn2d/`: images and settings for the synthetic tomography module

The data sets were provided by their authors for use with SEATREE.
Please cite the original publications, named in the subdirectories,
when using them.
