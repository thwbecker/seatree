# SEATREE packaged builds

Scripts to build a self-contained AppImage (Linux) or an application
bundle and DMG (macOS) of the SEATREE GUI, written by Dunyu Liu.
Both expect the compiled codes to have been built first with
`make all` in the repository root.

## Quick Build

From root directory:
```bash
./build.sh              # Build seatree.exe (main)
./build.sh clean        # Clean all build artifacts
```

## Build Scripts

- `build-appimage-gtk.sh` - Main AppImage with GTK4 dependencies
- `build.seatree.on.macos.sh` - macOS application bundle and DMG
- `clean.sh` - Remove all generated files

The portable AppImage script mentioned in earlier versions is not in
the repository at the moment. The AppImage scripts expect HC and
ConMan to have been built with `make all` in the root directory; they
have not been re-tested after the switch to git submodules.

## Output

- `seatree.exe` - Main AppImage (works on modern Linux)

## Requirements

Build tools are downloaded automatically:
- linuxdeploy-x86_64.AppImage
- appimagetool-x86_64.AppImage
- linuxdeploy-plugin-gtk.sh

## Usage

Built AppImages are self-contained and require no installation:
```bash
./seatree.exe
```