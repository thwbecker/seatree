#!/bin/bash

# SEATREE Build Cleanup Script
# Removes all generated files and build artifacts from the distribute folder

set -e

echo "SEATREE Build Cleanup"
echo "===================="

# Function to safely remove files/directories
safe_remove() {
    local item="$1"
    if [ -e "$item" ]; then
        echo "Removing: $item"
        rm -rf "$item"
    else
        echo "Not found: $item (already clean)"
    fi
}

echo ""
echo "Cleaning build artifacts..."

# Remove build directories
safe_remove "SEATREE.AppDir"
safe_remove "SEATREE-GTK.AppDir"
safe_remove "SEATREE-Portable.AppDir"

echo ""
echo "Cleaning output AppImages..."

# Remove output AppImages
safe_remove "seatree.exe"
safe_remove "seatree-gtk.exe"
safe_remove "seatree-portable.exe"
safe_remove "SEATREE-*.AppImage"

echo ""
echo "Cleaning downloaded tools (optional)..."

# Ask user if they want to remove downloaded tools
read -p "Remove downloaded build tools? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    safe_remove "appimagetool-x86_64.AppImage"
    safe_remove "linuxdeploy-x86_64.AppImage"
    safe_remove "linuxdeploy-plugin-gtk-x86_64.AppImage"
    echo "Downloaded tools removed (will be re-downloaded on next build)"
else
    echo "Downloaded tools kept (faster next build)"
fi

echo ""
echo "Cleaning temporary files..."

# Remove any temporary files
safe_remove "*.tmp"
safe_remove "*.log"

echo ""
echo "Cleanup complete!"
echo ""
echo "What was cleaned:"
echo "  ✓ Build directories (SEATREE*.AppDir)"
echo "  ✓ Output AppImages (*.exe, *.AppImage)"
echo "  ✓ Temporary files"
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "  ✓ Downloaded build tools"
else
    echo "  - Downloaded build tools (kept)"
fi
echo ""
echo "Ready for a fresh build!"