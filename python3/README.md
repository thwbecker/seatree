# SEATREE Python3 Application

GTK4-based GUI application for solid earth analysis and modeling.

## Structure

- `seatree/` - Main application package
  - `gui/` - GTK4 GUI components
  - `modules/` - Scientific computing modules (HC, ConMan, syn2d, larry, etc.)
  - `plotter/` - GMT and matplotlib plotting backends
  - `xml/` - Configuration file handling
- `conf/` - Module configuration files (XML); `conf/install/` lists the modules and their installers, the module `*Conf.xml` files hold paths relative to this directory

## Running

From the root directory:
```bash
./seatree              # or: python3 run.seatree.py
```

## Requirements

- Python 3.x with PyGObject (gi) for GTK 4, NumPy, SciPy, matplotlib (3.9 or newer)
- GMT 6 for plotting (GMT 4 is an optional legacy route)
- NetCDF libraries
- the compiled module codes, see the top-level README.md and Makefile

## Modules

- **HC** - Mantle circulation modeling
- **ConMan** - Convection modeling
- **syn2d** - 2D seismic tomography
- **larry/larry3d** - Block modeling
- **NonLinLoc** - Earthquake location