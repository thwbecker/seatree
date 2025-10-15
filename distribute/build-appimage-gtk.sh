#!/bin/bash

# Build script for SEATREE AppImage using linuxdeploy with GTK plugin
# This creates a truly self-contained AppImage with all GTK4 dependencies

set -e  # Exit on any error

echo "Building SEATREE self-contained AppImage with GTK4..."

# Function to check if required tools are available
check_tools() {
    echo "Checking required tools..."

    if [ ! -f "linuxdeploy-x86_64.AppImage" ]; then
        echo "linuxdeploy-x86_64.AppImage not found. Downloading..."
        wget -q --show-progress https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage
        if [ $? -ne 0 ]; then
            echo "Error: Failed to download linuxdeploy-x86_64.AppImage"
            exit 1
        fi
        chmod +x linuxdeploy-x86_64.AppImage
        echo "Downloaded and made executable: linuxdeploy-x86_64.AppImage"
    fi

    if [ ! -f "linuxdeploy-plugin-gtk.sh" ]; then
        echo "linuxdeploy-plugin-gtk.sh not found. Downloading..."
        wget -q --show-progress https://raw.githubusercontent.com/linuxdeploy/linuxdeploy-plugin-gtk/master/linuxdeploy-plugin-gtk.sh
        if [ $? -ne 0 ]; then
            echo "Error: Failed to download linuxdeploy-plugin-gtk.sh"
            exit 1
        fi
        chmod +x linuxdeploy-plugin-gtk.sh
        echo "Downloaded and made executable: linuxdeploy-plugin-gtk.sh"
    fi

    # Check if system has GTK4 installed
    if ! pkg-config --exists gtk4; then
        echo "Error: GTK4 development libraries not found on build system."
        echo "Please install: sudo apt-get install libgtk-4-dev"
        exit 1
    fi

    echo "All required tools found."
}

# Function to prepare application directory
prepare_appdir() {
    local appdir="$1"
    echo "Preparing AppDir structure: $appdir"

    # Clean previous build
    rm -rf "$appdir"

    # Create directory structure
    mkdir -p "$appdir/usr"/{bin,lib,share}
    mkdir -p "$appdir/usr/share"/{applications,icons/hicolor/256x256/apps}

    # Copy application files
    echo "Copying SEATREE application files..."
    mkdir -p "$appdir/usr/share/seatree"
    cp -r ../python3 "$appdir/usr/share/seatree/"
    cp -r ../modules "$appdir/usr/share/seatree/"

    # Bundle Python dependencies more comprehensively
    echo "Bundling Python dependencies..."

    # Create Python site-packages directory that matches system version
    PYTHON_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
    mkdir -p "$appdir/usr/lib/python$PYTHON_VERSION/site-packages"

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
                cp -r "$path" "$appdir/usr/lib/python$PYTHON_VERSION/site-packages/"
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

    # Also try to copy from multiple possible system locations
    # This ensures we get all the native extensions
    for sys_path in /usr/lib/python*/dist-packages /usr/local/lib/python*/dist-packages; do
        if [ -d "$sys_path/gi" ]; then
            echo "Also copying gi from: $sys_path"
            cp -r "$sys_path/gi"/* "$appdir/usr/lib/python$PYTHON_VERSION/site-packages/gi/" 2>/dev/null || true
        fi
    done

    # Copy GMT and NetCDF if they exist
    if [ -d "../gmt-4.5.18" ]; then
        echo "Copying GMT..."
        cp -r ../gmt-4.5.18 "$appdir/usr/share/seatree/"
    fi

    if [ -d "../netcdf-c-4.9.3-rc1" ]; then
        echo "Copying NetCDF..."
        cp -r ../netcdf-c-4.9.3-rc1 "$appdir/usr/share/seatree/"
    fi

    # Update configuration paths for AppImage
    echo "Updating configuration paths..."
    if [ -d "$appdir/usr/share/seatree/python3/conf" ]; then
        # Replace paths using suffix matching - works regardless of installation location
        # This matches paths like /any/path/seatree.dev/seatree.dev/... regardless of prefix
        for xmlfile in "$appdir/usr/share/seatree/python3/conf"/*.xml "$appdir/usr/share/seatree/python3/conf"/*/*.xml; do
            if [ -f "$xmlfile" ]; then
                sed -i 's|[^<>]*seatree\.dev/seatree\.dev/python3/seatree/modules/|/tmp/__APPIMAGE_RUNTIME__/python3/seatree/modules/|g; s|[^<>]*seatree\.dev/seatree\.dev/modules/|/tmp/__APPIMAGE_RUNTIME__/modules/|g; s|[^<>]*seatree\.dev/seatree\.dev/gmt-[^<>/]*/bin|/tmp/__APPIMAGE_RUNTIME__/gmt-4.5.18/bin|g' "$xmlfile"
            fi
        done

        echo "Replaced hardcoded paths with placeholders using suffix matching"
    fi

    # Create the netcdf symlinks during build
    if [ -f "$appdir/usr/share/seatree/netcdf-c-4.9.3-rc1/lib/libnetcdf.so.22" ]; then
        ln -sf libnetcdf.so.22 "$appdir/usr/share/seatree/netcdf-c-4.9.3-rc1/lib/libnetcdf.so.7"
        ln -sf libnetcdf.so.22 "$appdir/usr/share/seatree/netcdf-c-4.9.3-rc1/lib/libnetcdf.so.19"
    fi
}

# Function to create launcher script
create_launcher() {
    local appdir="$1"
    echo "Creating main launcher script..."

    cat > "$appdir/usr/bin/seatree-launcher" << 'EOF'
#!/bin/bash

HERE="$(dirname "$(dirname "$(readlink -f "${0}")")")"

# Set up environment variables for SEATREE
export APPDIR="$HERE"
export SEATREE="$HERE/usr/share/seatree"
export GMT4HOME="$HERE/usr/share/seatree/gmt-4.5.18"
export GMTHOME="$GMT4HOME"
export GMT_GSHHG_DATA="$GMT4HOME/gshhg-gmt-2.3.7"
export NETCDFHOME="$HERE/usr/share/seatree/netcdf-c-4.9.3-rc1"
export ARCH="x86_64"

# Set up GMT environment variables
export GMT_SHAREDIR="$GMT4HOME/share"
export GMT_DATADIR="$GMT4HOME/share"
export GMT_FONTPATH="$GMT4HOME/share/pslib"

# Set up library paths - NetCDF lib MUST come first for GMT binaries to find libnetcdf.so.19
export LD_LIBRARY_PATH="$NETCDFHOME/lib:$GMT4HOME/lib:$HERE/usr/lib:$LD_LIBRARY_PATH"
export PATH="$HERE/usr/bin:$GMT4HOME/bin:$NETCDFHOME/bin:$PATH"

# Set up Python path
PYTHON_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
export PYTHONPATH="$HERE/usr/lib/python$PYTHON_VERSION/site-packages:$HERE/usr/share/seatree/python3:$HERE/usr/share/seatree:$PYTHONPATH"

# Set APPIMAGE_ROOT for path resolution in config files
export APPIMAGE_ROOT="$HERE/usr/share/seatree"

# Replace placeholder paths in config files at runtime
if [ -d "$HERE/usr/share/seatree/python3/conf" ]; then
    # Replace the placeholder
    find "$HERE/usr/share/seatree/python3/conf" -name "*.xml" -exec sed -i \
        "s|/tmp/__APPIMAGE_RUNTIME__|$HERE/usr/share/seatree|g" {} \; 2>/dev/null || true
    # Replace ANY absolute path that ends with the expected directory structure
    # This catches paths regardless of where they were installed on build machine
    find "$HERE/usr/share/seatree/python3/conf" -name "*.xml" -exec sed -i \
        "s|[^<>]*seatree\\.dev/seatree\\.dev/python3/seatree/modules/|$HERE/usr/share/seatree/python3/seatree/modules/|g; \
         s|[^<>]*seatree\\.dev/seatree\\.dev/modules/|$HERE/usr/share/seatree/modules/|g; \
         s|[^<>]*seatree\\.dev/seatree\\.dev/gmt-[^<>/]*/bin|$HERE/usr/share/seatree/gmt-4.5.18/bin|g" {} \; 2>/dev/null || true
fi

# Ensure GTK can find its modules and themes
export GTK_PATH="$HERE/usr/lib/gtk-4.0"
export GTK_DATA_PREFIX="$HERE/usr"
export GTK_EXE_PREFIX="$HERE/usr"

# Ensure HOME is set
if [ -z "$HOME" ]; then
    export HOME="$USER_HOME"
fi
if [ -z "$HOME" ]; then
    export HOME="/tmp"
fi

# Use X11 backend to avoid Wayland issues
export GDK_BACKEND=x11

# Launch the application
cd "$HERE/usr/share/seatree"
exec python3 "$HERE/usr/share/seatree/python3/seatree/gui/SEATREE.py" "$@"
EOF

    chmod +x "$appdir/usr/bin/seatree-launcher"
}

# Function to create desktop file
create_desktop_file() {
    local appdir="$1"
    echo "Creating desktop file..."

    cat > "$appdir/seatree.desktop" << 'EOF'
[Desktop Entry]
Type=Application
Name=SEATREE
Comment=Solid Earth Analysis Tool for Research and Education
Exec=python3
Icon=seatree
StartupNotify=true
NoDisplay=false
Categories=Science;Education;
Keywords=geophysics;earth;science;modeling;
EOF
}

# Function to create application icon
create_icon() {
    local appdir="$1"
    echo "Creating application icon..."

    # Create a proper 256x256 PNG icon using ImageMagick if available
    if command -v convert >/dev/null 2>&1; then
        convert -size 256x256 xc:blue -fill white -gravity center -pointsize 72 -annotate +0+0 "S" "$appdir/seatree.png"
    else
        # Create a minimal valid 256x256 PNG icon
        python3 -c "
from PIL import Image
import os
img = Image.new('RGBA', (256, 256), (0, 100, 200, 255))
img.save('$appdir/seatree.png')
" 2>/dev/null || {
        # Fallback: create using dd (creates a valid but corrupt PNG)
        echo -e '\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x01\x00\x00\x00\x01\x00\x08\x02\x00\x00\x00\x90wS\xde\x00\x00\x00\x0cIDATx\x9cc\xf8\x0f\x00\x01\x01\x01\x00\x18\xdd\x8d\xb4\x00\x00\x00\x00IEND\xaeB`\x82' > "$appdir/seatree.png"
        }
    fi

    # Copy to standard icon location
    cp "$appdir/seatree.png" "$appdir/usr/share/icons/hicolor/256x256/apps/seatree.png" 2>/dev/null || true
}

# Function to create AppRun script
create_apprun() {
    local appdir="$1"
    echo "Creating AppRun script..."

    cat > "$appdir/AppRun" << 'EOF'
#!/bin/bash

HERE="$(dirname "$(readlink -f "${0}")")"

# Set up environment variables for SEATREE
export APPDIR="$HERE"
export SEATREE="$HERE/usr/share/seatree"
export GMT4HOME="$HERE/usr/share/seatree/gmt-4.5.18"
export GMTHOME="$GMT4HOME"
export GMT_GSHHG_DATA="$GMT4HOME/gshhg-gmt-2.3.7"
export NETCDFHOME="$HERE/usr/share/seatree/netcdf-c-4.9.3-rc1"
export ARCH="x86_64"

# Set up GMT environment variables
export GMT_SHAREDIR="$GMT4HOME/share"
export GMT_DATADIR="$GMT4HOME/share"
export GMT_FONTPATH="$GMT4HOME/share/pslib"

# Set up library paths - NetCDF lib MUST come first for GMT binaries to find libnetcdf.so.19
export LD_LIBRARY_PATH="$NETCDFHOME/lib:$GMT4HOME/lib:$HERE/usr/lib:$LD_LIBRARY_PATH"
export PATH="$HERE/usr/bin:$GMT4HOME/bin:$NETCDFHOME/bin:$PATH"

# Set up Python path
PYTHON_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
export PYTHONPATH="$HERE/usr/lib/python$PYTHON_VERSION/site-packages:$HERE/usr/share/seatree/python3:$HERE/usr/share/seatree:$PYTHONPATH"

# Set APPIMAGE_ROOT for path resolution in config files
export APPIMAGE_ROOT="$HERE/usr/share/seatree"

# Replace placeholder paths in config files at runtime
if [ -d "$HERE/usr/share/seatree/python3/conf" ]; then
    # Replace the placeholder
    find "$HERE/usr/share/seatree/python3/conf" -name "*.xml" -exec sed -i \
        "s|/tmp/__APPIMAGE_RUNTIME__|$HERE/usr/share/seatree|g" {} \; 2>/dev/null || true
    # Replace ANY absolute path that ends with the expected directory structure
    # This catches paths regardless of where they were installed on build machine
    find "$HERE/usr/share/seatree/python3/conf" -name "*.xml" -exec sed -i \
        "s|[^<>]*seatree\\.dev/seatree\\.dev/python3/seatree/modules/|$HERE/usr/share/seatree/python3/seatree/modules/|g; \
         s|[^<>]*seatree\\.dev/seatree\\.dev/modules/|$HERE/usr/share/seatree/modules/|g; \
         s|[^<>]*seatree\\.dev/seatree\\.dev/gmt-[^<>/]*/bin|$HERE/usr/share/seatree/gmt-4.5.18/bin|g" {} \; 2>/dev/null || true
fi

# Ensure GTK can find its modules and themes
export GTK_PATH="$HERE/usr/lib/gtk-4.0"
export GTK_DATA_PREFIX="$HERE/usr"
export GTK_EXE_PREFIX="$HERE/usr"

# Ensure HOME is set
if [ -z "$HOME" ]; then
    export HOME="$USER_HOME"
fi
if [ -z "$HOME" ]; then
    export HOME="/tmp"
fi

# Use X11 backend to avoid Wayland issues
export GDK_BACKEND=x11

# Launch the application
cd "$HERE/usr/share/seatree"
exec "$HERE/usr/bin/python3" "$HERE/usr/share/seatree/python3/seatree/gui/SEATREE.py" "$@"
EOF

    chmod +x "$appdir/AppRun"
}

# Main build function
main() {
    local appdir="SEATREE-GTK.AppDir"
    local output_name="seatree.exe"

    # Ensure we use system Python, not conda
    echo "Deactivating conda and using system Python..."
    unset CONDA_DEFAULT_ENV
    unset CONDA_PREFIX
    export PATH="/usr/bin:/bin:/usr/sbin:/sbin:${PATH}"

    echo "Using Python: $(which python3)"
    echo "Python version: $(python3 --version)"

    # Check prerequisites
    check_tools

    # Prepare application directory
    prepare_appdir "$appdir"

    # Create launcher and desktop files
    create_launcher "$appdir"
    create_desktop_file "$appdir"
    create_icon "$appdir"
    create_apprun "$appdir"

    # Set environment variables for linuxdeploy
    export LINUXDEPLOY="$(pwd)/linuxdeploy-x86_64.AppImage"
    export DEPLOY_GTK_VERSION=4
    export GDK_BACKEND=x11
    export ARCH=x86_64

    echo "Running linuxdeploy with GTK plugin..."

    # Download and bundle a portable Python3 if needed
    echo "Setting up portable Python3..."

    # Create a wrapper script that uses bundled Python libs
    cat > "$appdir/usr/bin/python3-portable" << 'PYTHON_EOF'
#!/bin/bash
HERE="$(dirname "$(dirname "$(readlink -f "${0}")")")"
export LD_LIBRARY_PATH="$HERE/usr/lib:$LD_LIBRARY_PATH"
exec "$HERE/usr/bin/python3" "$@"
PYTHON_EOF
    chmod +x "$appdir/usr/bin/python3-portable"

    # Use linuxdeploy with GTK plugin to bundle GTK4 dependencies
    # Use the system python3 but with library bundling
    ./linuxdeploy-x86_64.AppImage \
        --appdir "$appdir" \
        --executable "/usr/bin/python3" \
        --desktop-file "$appdir/seatree.desktop" \
        --icon-file "$appdir/seatree.png" \
        --plugin gtk \
        --output appimage

    # The output will be named based on the desktop file, let's rename it
    appimage_file=$(ls SEATREE-*.AppImage 2>/dev/null | head -1)
    if [ -n "$appimage_file" ] && [ -f "$appimage_file" ]; then
        mv "$appimage_file" "$output_name"
        echo "Self-contained AppImage built successfully: $output_name"
        echo "This AppImage includes all GTK4 dependencies and should run on any Linux system!"
        echo "File size: $(du -h "$output_name" | cut -f1)"
    else
        echo "Error: AppImage was not created successfully"
        exit 1
    fi
}

# Run the main function
main "$@"