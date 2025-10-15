#!/bin/bash

# SEATREE AppImage Build Script
# This script provides easy access to build AppImages from the root directory

set -e

echo "SEATREE AppImage Builder"
echo "======================="

case "${1:-main}" in
    "main"|"")
        echo "Building main AppImage (seatree.exe)..."
        cd distribute && ./build-appimage-gtk.sh
        ;;
    "portable")
        echo "Building portable AppImage (seatree-portable.exe)..."
        cd distribute && ./build-appimage-gtk-portable.sh
        ;;
    "clean")
        echo "Cleaning build artifacts..."
        cd distribute && ./clean.sh
        ;;
    "help"|"-h"|"--help")
        echo "Usage: $0 [main|portable|clean|help]"
        echo ""
        echo "Options:"
        echo "  main      - Build main AppImage (seatree.exe) [default]"
        echo "  portable  - Build portable AppImage (seatree-portable.exe)"
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