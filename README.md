# SEATREE - Solid Earth Teaching and Research Environment

**Thorsten Becker** - thwbecker@post.harvard.edu - 10/2024

SEATREE is a modular and user-friendly software to facilitate using solid Earth research tools in the classroom and for interdisciplinary scientific collaboration.

We use python wrappers and make use of modern software design concepts, while remaining compatible with traditional scientific coding. Our goals are to provide a fully contained, yet transparent package that lets users operate in an easy, graphically supported "black box" mode, while also allowing to look under the hood. In the long run, we envision SEATREE to contribute to new ways of sharing scientific research, and making (numerical) experiments truly reproducible again.

## Modules

SEATREE is module based, and the current version includes tools for:

- **2-D Mantle (Thermal) Convection** - Computing thermal convection models
- **3-D Body Wave Mantle Seismic Tomography** - Seismic imaging of the mantle
- **3-D Spherical Mantle Flow** - A driver for the HC tool
- **Surface Wave Phase Velocity Tomography** - Inverting for Earth structure
- **2-D Synthetic Tomography** - Teaching module
- **Earthquake Location Inversions** - Rudimentary module for earthquake location

The main software design consists of transparent python wrappers that drive the modules, including a GMT plotting tool, a VTK/Paraview 3-D visualization interface, and a graphical user interface.

## Installation

### Quick Start (Ubuntu 22)

```bash
# Clone SEATREE from GitHub repository
git clone https://github.com/thwbecker/seatree/
cd seatree

# If Ubuntu packages are not installed, use -e option to install packages, netcdf, gmt, and seatree
./install.seatree.sh -e ubuntu

# If Ubuntu packages are already installed, use -m option to install netcdf, gmt4, and seatree
./install.seatree.sh -m ubuntu
```

### Running SEATREE

```bash
python3 run.seatree.py  # to start the SEATREE APP
```

### Other Installation Options

For VirtualBox installation and other download options, see:

http://www-udc.ig.utexas.edu/external/becker/seatree/wiki/download.html

## Resources

**Project Website:** http://www-udc.ig.utexas.edu/external/becker/seatree/


## Reference

Milner, K., Becker, T. W., Boschi, L., Sain, J., Schorlemmer, D. and H. Waterhouse: The Solid Earth Research and Teaching Environment: a new software framework to share research tools in the classroom and across disciplines. *Eos Trans. AGU*, 90, 12, 2009.

http://www-udc.ig.utexas.edu/external/becker/preprints/mbbssw08.pdf

## License

SEATREE is freely available under the GNU General Public License v2. See [LICENSE.md](LICENSE.md) for the full license text.

## Support and Contributing

If you want to use SEATREE in a classroom setting, we might be able to offer you some installation support and always welcome your feedback. Also, if you like to add your own module to SEATREE, please let us know; we might be able to provide some assistance.

## Acknowledgments

Partial support was provided through NSF-CAREER.

## Third-Party Components and Copyright

SEATREE includes and may be distributed with the following third-party components:

- **fstrack**: Set of routines to trace particles in mantle flow. Routines use precomputed velocities on GMT grd files. See: Becker et al. (GJI, 155, 696, 2003).

- **d-rex**: DREX Kaminski & Ribe texture routines by Eduard Kaminski and Jules Browaeys. See: Kaminski et al. (GJI, 157, 1, 2004) and Browaeys & Chevrot (GJI, 159, 667, 2004).

- **single_layer and multi_layer**: Shear wave splitting routines by Vera Schulte-Pelkum. See: Schulte-Pelkum & Blackman (GJI, 154, 166, 2003).

- **eispack**: EISPACK eigensystem routines

- **menke_splitting**: Bill Menke's cross-correlation routines

Copyright for these components remains with the original authors.

For the rest: Copyright (c) Thorsten Becker 2004-2005
