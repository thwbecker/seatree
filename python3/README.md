# SEATREE python package

The GTK 4 application and the module drivers. Start it from the
repository root with `./seatree` (or `python3 run.seatree.py`), which
sets the environment variables the modules expect and runs
`seatree/gui/SEATREE.py`.

## Layout

- `seatree/gui/`: application and main window (`SEATREE.py`), start
  dialog, save/load dialogs, the script viewer
- `seatree/modules/`: one package per research code; `module.py` is
  the base class, `example/` a minimal template
- `seatree/plotter/`: `gmt/` (PostScript through GMT, converted to
  PNG for display) and `matPlotLib/` (drawn into the window) plotters,
  plus the image widget they share
- `seatree/gmt/`: wrapper that builds GMT command lines
- `seatree/util/`: subprocess runner (`scriptRunner.py`), PostScript
  conversion, small helpers
- `seatree/xml/`: reading and writing of the configuration files
- `seatree/install/`: the installer (`install.py`) and helpers
- `conf/`: configuration; see below
- `data/`: input data offered by the modules
- `img/`: images used by the GUI

## Module interface

A module is a class derived from `seatree.modules.module.Module`.
The main window calls:

- `setDefaults(mainWindow)`: receives the main window, from which the
  module takes the temporary file prefix, the GMT path and the
  PostScript converter, and creates its plotter
- `getPanel(mainWindow)`: returns the GTK widget with the module's
  controls (by convention a `<Name>GUI` class in its own file)
- `getPlotter()`: returns the plotter whose widget is shown on the
  right
- `getSettings()` and `loadSettings(node)`: the module's state as an
  XML element, for Save/Load Settings
- `updatePlot()`: redraw after a plotter setting changed
- `cleanup()`: called when the module is unloaded or SEATREE exits

Modules run the compiled codes as subprocesses with
`ScriptRunner`, in a temporary directory under `/tmp/seatree.*`, and
append every command to their output log, which View Script shows.
The commands reference the binaries by the paths in the module's
configuration file.

## Installer

`seatree/install/install.py` (run by `../install/configure.python3.gtk4`)
reads `conf/install/*.xml` in order, imports the installer class named
there for each module (`installImportName`, `installClassName`),
and calls its `install()`. The installer classes build the module's
code through the top-level Makefile (`make -C <root> <target>`,
`installUtils.buildWithMake`) and write the module configuration
file. `install.py` then writes `conf/conf.xml` and checks for GMT and
ImageMagick.

## Configuration files

- `conf/conf.xml`: the list of modules (import name, class name,
  directory) and the GMT and `convert` paths. An empty `gmtPath`
  means the GMT in the system path; a configured path that does not
  exist falls back to the system path with a message.
- `conf/<module>/<Module>Conf.xml`: where each module finds its
  binaries, for example `hcPath`, `conmanPath` (the `src/` directory,
  the binary is one level up), `larryPath`, `binPath` for NonLinLoc.
- `conf/install/<n>_<module>.xml`: the installer registry.
- `conf/baseUI.xml`: the menu definition.

Relative paths in these files are relative to this `python3/`
directory. The committed defaults match an unmoved checkout after
`make all`; the installer rewrites them if asked.

## Adding a module

1. Create `seatree/modules/<name>/` with the driver class, the GUI
   class, and an installer class (start from `example/`).
2. Add `conf/install/<n>_<name>.xml` with the import and class names;
   the number sets the order in the start dialog.
3. Put the compiled code under `../modules/` and add a target for it
   to the top-level Makefile; have the installer call
   `installUtils.buildWithMake(pythonRoot, "<target>")`.
4. Run `../install/configure.python3.gtk4` to regenerate `conf/conf.xml`.

## Requirements

Python 3 with PyGObject for GTK 4, NumPy, SciPy and matplotlib 3.9 or
newer; GMT 6, Ghostscript and ImageMagick for the GMT based modules.
See the top-level README for packages and the build of the compiled
codes.
