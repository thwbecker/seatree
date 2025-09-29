# SEATREE Python3 Application

GTK4-based GUI application for solid earth analysis and modeling.

## Structure

- `seatree/` - Main application package
  - `gui/` - GTK4 GUI components
  - `modules/` - Scientific computing modules (HC, ConMan, syn2d, larry, etc.)
  - `plotter/` - GMT and matplotlib plotting backends
  - `xml/` - Configuration file handling
- `conf/` - Module configuration files (XML)

## Running

From root directory:
```bash
python3 run.seatree.python3.gtk4.py
```

## Requirements

- Python 3.x with PyGObject (gi)
- GTK4 development libraries
- GMT 4.x for plotting
- NetCDF libraries

## Modules

- **HC** - Mantle circulation modeling
- **ConMan** - Convection modeling
- **syn2d** - 2D seismic tomography
- **larry/larry3d** - Block modeling
- **NonLinLoc** - Earthquake location