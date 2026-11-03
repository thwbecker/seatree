#!/bin/bash

# SEATREE Distribution Build Script
# Builds platform-specific distributions (AppImage for Linux, .app/.dmg for macOS)

set -e

# Detect platform
PLATFORM="$(uname -s)"

if [[ "$PLATFORM" == "Darwin" ]]; then
    echo "SEATREE macOS Builder"
    echo "====================="

    case "${1:-main}" in
        "main"|"")
            echo "Building macOS app bundle and DMG..."
            cd distribute && ./build.seatree.on.macos.sh
            ;;
        "clean")
            echo "Cleaning build artifacts..."
            rm -rf distribute/dist
            ;;
        "help"|"-h"|"--help")
            echo "Usage: $0 [main|clean|help]"
            echo ""
            echo "Options:"
            echo "  main   - Build macOS .app and .dmg [default]"
            echo "  clean  - Clean all build artifacts"
            echo "  help   - Show this help message"
            echo ""
            echo "Output: distribute/dist/Seatree.app and Seatree-YYYYMMDD.dmg"
            ;;
        *)
            echo "Error: Unknown option '$1'"
            echo "Use '$0 help' for usage information"
            exit 1
            ;;
    esac
elif [[ "$PLATFORM" == "Linux" ]]; then
    echo "SEATREE AppImage Builder"
    echo "======================="

    case "${1:-main}" in
        "main"|"")
            echo "Building main AppImage (seatree.exe)..."
            cd distribute && ./build-appimage-gtk.sh
            ;;
        "portable")
            echo "Error: the portable AppImage build script (build-appimage-gtk-portable.sh)"
            echo "is not part of this repository at the moment."
            exit 1
            ;;
        "clean")
            echo "Cleaning build artifacts..."
            cd distribute && ./clean.sh
            ;;
        "help"|"-h"|"--help")
            echo "Usage: $0 [main|clean|help]"
            echo ""
            echo "Options:"
            echo "  main      - Build main AppImage (seatree.exe) [default]"
            echo "  clean     - Clean all build artifacts and temporary files"
            echo "  help      - Show this help message"
            echo ""
            echo "Build files are located in the 'distribute/' folder"
            ;;
        *)
            echo "Error: Unknown option '$1'"
            echo "Use '$0 help' for usage information"
            exit 1
            ;;
    esac
else
    echo "ERROR: Unsupported platform: $PLATFORM"
    exit 1
fi