# SEATREE AppImage Distribution

Build self-contained AppImage executables for SEATREE GUI application.

## Quick Build

From root directory:
```bash
./build.sh              # Build seatree.exe (main)
./build.sh portable     # Build seatree-portable.exe (older systems)
./build.sh clean        # Clean all build artifacts
```

## Build Scripts

- `build-appimage-gtk.sh` - Main AppImage with GTK4 dependencies
- `build-appimage-gtk-portable.sh` - Portable version for older Linux systems
- `clean.sh` - Remove all generated files

## Output

- `seatree.exe` - Main AppImage (works on modern Linux)
- `seatree-portable.exe` - Portable AppImage (works on older GLIBC)

## Requirements

Build tools are downloaded automatically:
- linuxdeploy-x86_64.AppImage
- appimagetool-x86_64.AppImage
- linuxdeploy-plugin-gtk.sh

## Usage

Built AppImages are self-contained and require no installation:
```bash
./seatree.exe
./seatree-portable.exe
```