#!/bin/bash

# Build SEATREE AppImage with portable Python for GLIBC compatibility
# This version works on older Linux systems

set -e

echo "Building SEATREE AppImage with portable Python..."

# Clean and create AppDir
rm -rf SEATREE-Portable.AppDir
mkdir -p SEATREE-Portable.AppDir/usr/{bin,lib,share/seatree}

# Use system Python but bundle it properly
echo "Using system Python for portable build..."

# Create Python site-packages directory that matches system version
PYTHON_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
mkdir -p "SEATREE-Portable.AppDir/usr/lib/python$PYTHON_VERSION/site-packages"

# Function to copy Python package with all dependencies
copy_python_package() {
    local pkg_name="$1"
    local source_paths=(
        "/usr/lib/python3/dist-packages/$pkg_name"
        "/usr/lib/python$PYTHON_VERSION/dist-packages/$pkg_name"
        "/usr/local/lib/python3/dist-packages/$pkg_name"
        "/usr/local/lib/python$PYTHON_VERSION/dist-packages/$pkg_name"
    )

    for path in "${source_paths[@]}"; do
        if [ -d "$path" ]; then
            echo "Found $pkg_name at: $path"
            cp -r "$path" "SEATREE-Portable.AppDir/usr/lib/python$PYTHON_VERSION/site-packages/"
            return 0
        fi
    done
    echo "Warning: $pkg_name not found in standard locations"
    return 1
}

# Copy essential packages
copy_python_package "gi"
copy_python_package "cairo"
copy_python_package "PIL" || copy_python_package "Pillow"

# Also copy from system locations to ensure we get all native extensions
for sys_path in /usr/lib/python*/dist-packages /usr/local/lib/python*/dist-packages; do
    if [ -d "$sys_path/gi" ]; then
        echo "Also copying gi from: $sys_path"
        cp -r "$sys_path/gi"/* "SEATREE-Portable.AppDir/usr/lib/python$PYTHON_VERSION/site-packages/gi/" 2>/dev/null || true
    fi
done

echo "Python packages copied successfully"

# Copy application files
echo "Copying SEATREE application files..."
cp -r ../python3 SEATREE-Portable.AppDir/usr/share/seatree/
cp -r ../modules SEATREE-Portable.AppDir/usr/share/seatree/

# Copy GMT and NetCDF
if [ -d "../gmt-4.5.18" ]; then
    cp -r ../gmt-4.5.18 SEATREE-Portable.AppDir/usr/share/seatree/
fi
if [ -d "../netcdf-c-4.9.3-rc1" ]; then
    cp -r ../netcdf-c-4.9.3-rc1 SEATREE-Portable.AppDir/usr/share/seatree/

    # Create NetCDF symlinks for version compatibility
    if [ -f "SEATREE-Portable.AppDir/usr/share/seatree/netcdf-c-4.9.3-rc1/lib/libnetcdf.so.22" ]; then
        ln -sf libnetcdf.so.22 SEATREE-Portable.AppDir/usr/share/seatree/netcdf-c-4.9.3-rc1/lib/libnetcdf.so.7
        ln -sf libnetcdf.so.22 SEATREE-Portable.AppDir/usr/share/seatree/netcdf-c-4.9.3-rc1/lib/libnetcdf.so.19
        echo "Created NetCDF version symlinks"
    fi
fi

# Update configuration paths
echo "Updating configuration paths..."
if [ -d "SEATREE-Portable.AppDir/usr/share/seatree/python3/conf" ]; then
    find SEATREE-Portable.AppDir/usr/share/seatree/python3/conf -name "*.xml" -exec sed -i 's|/home/utig5/dliu/seatree.dev/seatree.dev|APPIMAGE_ROOT|g' {} \\;
    find SEATREE-Portable.AppDir/usr/share/seatree/python3/conf -name "*.xml" -exec sed -i 's|/home/staff/dliu/seatree.dev/seatree.dev|APPIMAGE_ROOT|g' {} \\;
fi

# Create AppRun script
cat > SEATREE-Portable.AppDir/AppRun << 'EOF'
#!/bin/bash

HERE="$(dirname "$(readlink -f "${0}")")"

# Set up environment variables for SEATREE
export APPDIR="$HERE"
export SEATREE="$HERE/usr/share/seatree"
export GMT4HOME="$HERE/usr/share/seatree/gmt-4.5.18"
export GMTHOME="$GMT4HOME"
export GMT_GSHHG_DATA="$GMT4HOME/gshhg-gmt-2.3.7"
export NETCDFHOME="$HERE/usr/share/seatree/netcdf-c-4.9.3-rc1"

# Set up GMT environment variables
export GMT_SHAREDIR="$GMT4HOME/share"
export GMT_DATADIR="$GMT4HOME/share"
export GMT_FONTPATH="$GMT4HOME/share/pslib"

# Set up library paths
export LD_LIBRARY_PATH="$HERE/usr/lib:$NETCDFHOME/lib:$LD_LIBRARY_PATH"
export PATH="$HERE/usr/bin:$NETCDFHOME/lib:$GMT_GSHHG_DATA:$PATH"

# Set up Python path
export PYTHONPATH="$HERE/usr/share/seatree/python3:$HERE/usr/share/seatree:$PYTHONPATH"

# Set APPIMAGE_ROOT for path resolution
export APPIMAGE_ROOT="$HERE/usr/share/seatree"

# Use X11 backend
export GDK_BACKEND=x11

# Ensure HOME is set
if [ -z "$HOME" ]; then
    export HOME="/tmp"
fi

# Launch with system Python
cd "$HERE/usr/share/seatree"
exec python3 "$HERE/usr/share/seatree/python3/seatree/gui/SEATREE.py" "$@"
EOF

chmod +x SEATREE-Portable.AppDir/AppRun

# Create desktop file
cat > SEATREE-Portable.AppDir/seatree.desktop << 'EOF'
[Desktop Entry]
Type=Application
Name=SEATREE
Comment=Solid Earth Analysis Tool for Research and Education
Exec=AppRun
Icon=seatree
Categories=Science;Education;
EOF

# Create simple icon
echo "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8/5+hHgAHggJ/PchI7wAAAABJRU5ErkJggg==" | base64 -d > SEATREE-Portable.AppDir/seatree.png 2>/dev/null || echo "Could not create icon"

# Package as AppImage
echo "Creating portable AppImage..."
if [ ! -f "appimagetool-x86_64.AppImage" ]; then
    wget "https://github.com/AppImage/AppImageKit/releases/download/continuous/appimagetool-x86_64.AppImage"
    chmod +x appimagetool-x86_64.AppImage
fi

ARCH=x86_64 ./appimagetool-x86_64.AppImage SEATREE-Portable.AppDir seatree-portable.exe

echo "Portable AppImage created: seatree-portable.exe"
echo "This version should work on older Linux systems with older GLIBC!"