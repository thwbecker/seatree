#! /bin/bash

# The shell script is to set up environment for SEATREE.
# Currently supported systems include:
# 1. ubuntu: Ubuntu 22.04
# 2. docker: please see Dockerfile.
# 3. macos: MacOS 14

# GMT Version Selection (4 or 6)
# Default: GMT6 (system-installed via package manager)
# Optional: GMT4 (set GMTVERSION=4 to use legacy GMT4)
# Note: SEATREE now uses GMT6 by default. GMT4 is kept available
#       for users who need it, but GMT6 is recommended for new installations.
GMTVERSION=${GMTVERSION:-6}  # Default to 6 if not set

# Set up logging
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOGDIR="tmp"
mkdir -p "$LOGDIR"
LOGFILE="$LOGDIR/install_seatree_${TIMESTAMP}.log"

# Function to log messages to both screen and file
log_info() {
    echo "$@" | tee -a "$LOGFILE"
}

# Function to log only to file (suppress screen output)
log_verbose() {
    echo "$@" >> "$LOGFILE"
}

# Detect whether NetCDF is already available either from the bundled
# tree or from a system installation exposed via nc-config.
netcdf_is_installed() {
    if [ -n "$NETCDFHOME" ] && [ -d "$NETCDFHOME" ]; then
        if [ "$(uname -s)" = "Darwin" ]; then
            local c_home="${NETCDF_C_HOME:-$NETCDFHOME}"
            if [ -f "$c_home/lib/libnetcdf.dylib" ] && \
               [ -x "$c_home/bin/nc-config" ]; then
                return 0
            fi
        else
            if [ -f "$NETCDFHOME/lib/libnetcdf.so" ] && \
               [ -x "$NETCDFHOME/bin/nc-config" ]; then
                return 0
            fi
        fi
    fi

    if command -v nc-config >/dev/null 2>&1; then
        if nc-config --version >/dev/null 2>&1; then
            return 0
        fi
    fi

    return 1
}

detect_system_netcdf() {
    if command -v nc-config >/dev/null 2>&1; then
        local prefix
        prefix=$(nc-config --prefix 2>/dev/null | head -n 1)
        if [ -n "$prefix" ] && [ -d "$prefix" ]; then
            echo "$prefix"
            return 0
        fi
    fi
    return 1
}

log_info "=========================================="
log_info "SEATREE Installation Log"
log_info "Started: $(date)"
log_info "Log file: $LOGFILE"
log_info "GMT Version: ${GMTVERSION} (default: 6)"
log_info "=========================================="
log_info ""

while getopts "m:e:c:h" OPTION; do
    case $OPTION in
        m)
            MACH=$OPTARG
            ;;
        e)
            MACH=$OPTARG
            ENV="True"
            ;;
        h)
            echo "Usage: ./install.seatree.sh [-h] [-m Machine_name] [-e Machine_name] "
            echo "                                                                     "
            echo "Examples:                                                            "
            echo "                                                                     "
            echo "./install.seatree.sh -h                                 "
            echo " -----Display this help message                                      "
            echo "                                                                     "
            echo "./install.seatree.sh -m ubuntu                          "
            echo " -----Install SEATREE on ubuntu 22.04 (uses system GMT6 by default)  "
            echo "                                                                     "
            echo "GMTVERSION=4 ./install.seatree.sh -m ubuntu             "
            echo " -----Install SEATREE with GMT4 (optional, for legacy use cases)     "
            echo "                                                                     "
            echo "./install.seatree.sh -e macos                           "
            echo " -----Simply set up envs for SEATREE without installation            "
            echo " -----on macos                                                       "
            echo "                                                                     "
            echo "source ./install.seatree.sh                             "
            echo " -----Activate ENV VAR and add exes to PATH                          "
            echo "                                                                     "
            echo "Currently supported machines include:                                "
            echo " ubuntu/macos/docker                                                 "
            echo "                                                                     "
            echo "GMT Version:                                                         "
            echo " Default: GMT6 (system-installed via apt/homebrew)                   "
            echo " Optional: GMT4 (set GMTVERSION=4 before running this script)        "
            ;;
    esac
done 

if [ -n "$MACH" ]; then
    export MACHINE=$MACH
    if [ $MACHINE == "ubuntu" ]; then
        if [ -n "$ENV" ]; then
          sudo apt-get update
          sudo apt-get install git wget cmake python3 python3-pip python3-numpy
          sudo apt-get install python3-gi python3-gi-cairo gir1.2-gtk-4.0 libgtk-4-dev
          sudo apt-get install build-essential gfortran
          sudo apt-get install x11-apps libx11-dev libxt-dev libxaw7-dev ghostscript libhdf5-dev gawk
          # Install GMT6 (default) - required for SEATREE plotting modules
          sudo apt-get install gmt gmt-dcw gmt-gshhg libgmt-dev
          pip3 install --user matplotlib==3.9.2 #--break-system-packages # needs newer version of matplotlib to work.
        fi
    elif [ $MACHINE == "macos" ]; then
        echo "Installing SEATREE on MacOS 14 ... ..."
        if [ -n "$ENV" ]; then
            brew install gtk4 ghostscript gawk
            # Install GMT6 (default) - required for SEATREE plotting modules
            brew install gmt
            pip3 install --break-system-packages matplotlib==3.9.2 pygobject
        fi
    fi 

    export SEATREE=$(pwd)
    if [ "$GMTVERSION" == "4" ]; then
        export GMT4HOME=$(pwd)/gmt-4.5.18
        export GMTHOME=$GMT4HOME
    elif [ "$GMTVERSION" == "6" ]; then
        # For GMT6, detect installation location
        # First, try to find GMT6 using 'which gmt'
        if command -v gmt >/dev/null 2>&1; then
            GMT_BIN=$(which gmt)
            # Get the bin directory, then go up one level to get GMT home
            GMT_BIN_DIR=$(dirname "$GMT_BIN")
            export GMTHOME=$(dirname "$GMT_BIN_DIR")
            export GMT4HOME="$GMTHOME"  # For HC backward compatibility
            log_info "  -> Found GMT6 at: $GMT_BIN (GMTHOME=$GMTHOME)"
        elif [ -d "/usr/include/gmt" ] && [ -f "/usr/lib/x86_64-linux-gnu/libgmt.so" ]; then
            # Ubuntu/Debian system GMT6
            export GMTHOME="/usr"
            export GMT4HOME="$GMTHOME"  # For HC backward compatibility
        elif [ -d "/usr/local/include/gmt" ] && [ -f "/usr/local/lib/libgmt.dylib" ]; then
            # macOS Homebrew GMT6
            export GMTHOME="/usr/local"
            export GMT4HOME="$GMTHOME"
        else
            echo "ERROR: GMT6 not found in PATH, /usr, or /usr/local"
            echo "Please install GMT6 with: sudo apt-get install gmt gmt-dcw gmt-gshhg libgmt-dev"
            exit 1
        fi
    else
        echo "ERROR: GMTVERSION must be 4 or 6, got: $GMTVERSION"
        exit 1
    fi
    if [ -z "$NETCDFHOME" ]; then
        if [ "$GMTVERSION" = "6" ] && [ "$(uname -s)" = "Darwin" ] && command -v brew >/dev/null 2>&1 && brew ls --versions netcdf-fortran >/dev/null 2>&1; then
            export NETCDFHOME="$(brew --prefix netcdf-fortran)"
            log_info "Detected Homebrew netcdf-fortran: $NETCDFHOME"
            if brew ls --versions netcdf >/dev/null 2>&1; then
                export NETCDF_C_HOME="$(brew --prefix netcdf)"
                log_info "Detected Homebrew netcdf (C library): $NETCDF_C_HOME"
            fi
        else
            prefix=$(detect_system_netcdf) || true
            if [ -n "$prefix" ] && [ -d "$prefix" ]; then
                export NETCDFHOME="$prefix"
                log_info "Detected system NetCDF via nc-config: $NETCDFHOME"
            else
                export NETCDFHOME=$(pwd)/netcdf-c-4.9.3-rc1
                log_info "Using bundled NetCDF at $NETCDFHOME"
            fi
        fi
    else
        log_info "NETCDFHOME preset to $NETCDFHOME (honoring user override)"
    fi
    export ARCH=$(uname -m)

    log_info "[STEP 1/5] $(date +"%Y-%m-%d %H:%M:%S") - Checking NetCDF installation..."
    if netcdf_is_installed; then
        log_info "  -> NetCDF already installed (nc-config detected)."
    else
        log_info "  -> Installing netcdf-c-4.9.3-rc1..."
        bash install/install.netcdf.ubuntu22.sh >> "$LOGFILE" 2>&1
        if [ $? -eq 0 ]; then
            log_info "  -> NetCDF installation completed successfully."
        else
            log_info "  -> ERROR: NetCDF installation failed. Check $LOGFILE for details."
        fi
    fi

    log_info "[STEP 2/5] $(date +"%Y-%m-%d %H:%M:%S") - Checking GMT installation..."
    if [ "$GMTVERSION" == "4" ]; then
        if [ -e "gmt-4.5.18" ]; then
            log_info "  -> GMT 4 already installed, skipping."
        else
            log_info "  -> Installing GMT 4..."
            bash install/install.gmt4.ubuntu22.sh >> "$LOGFILE" 2>&1
            if [ $? -eq 0 ]; then
                log_info "  -> GMT 4 installation completed successfully."
            else
                log_info "  -> ERROR: GMT 4 installation failed. Check $LOGFILE for details."
            fi
        fi
    elif [ "$GMTVERSION" == "6" ]; then
        log_info "  -> Using system GMT 6 (skipping installation)"
    fi

    log_info "[STEP 3/5] $(date +"%Y-%m-%d %H:%M:%S") - Configuring Python3/GTK4..."
    yes '' | ./install/configure.python3.gtk4 >> "$LOGFILE" 2>&1
    if [ $? -eq 0 ]; then
        log_info "  -> Configuration completed successfully."
    else
        log_info "  -> ERROR: Configuration failed. Check $LOGFILE for details."
    fi

    log_info "[STEP 4/5] $(date +"%Y-%m-%d %H:%M:%S") - Installing ConMan v3.0.0..."
    bash install/install.conman.sh >> "$LOGFILE" 2>&1
    if [ $? -eq 0 ]; then
        log_info "  -> ConMan installation completed successfully."
    else
        log_info "  -> ERROR: ConMan installation failed. Check $LOGFILE for details."
    fi

    log_info "[STEP 5/5] $(date +"%Y-%m-%d %H:%M:%S") - Installing EQdyna..."
    if [ -n "$MACHINE" ]; then
        EQDYNA_MACHINE="$MACHINE"
    else
        system_name=$(uname | tr '[:upper:]' '[:lower:]')
        case "$system_name" in
            darwin*)
                EQDYNA_MACHINE="macos"
                ;;
            linux*)
                EQDYNA_MACHINE="ubuntu"
                ;;
            *)
                EQDYNA_MACHINE="$system_name"
                ;;
        esac
    fi
    bash modules/seismo/EQdyna/install-eqdyna.sh -m "$EQDYNA_MACHINE" >> "$LOGFILE" 2>&1
    if [ $? -eq 0 ]; then
        log_info "  -> EQdyna installation completed successfully."
    else
        log_info "  -> ERROR: EQdyna installation failed. Check $LOGFILE for details."
    fi

    log_info ""
    log_info "=========================================="
    log_info "Installation completed: $(date)"
    log_info "=========================================="
    log_info ""
    log_info "Running verification script..."
    python3 verify.seatree.py
fi

export SEATREE=$(pwd)
if [ "$GMTVERSION" == "4" ]; then
    export GMT4HOME=$(pwd)/gmt-4.5.18
    export GMTHOME=$GMT4HOME
    export GMT_GSHHG_DATA=$GMT4HOME/gshhg-gmt-2.3.7
elif [ "$GMTVERSION" == "6" ]; then
    # For GMT6, detect installation location
    # First, try to find GMT6 using 'which gmt'
    if command -v gmt >/dev/null 2>&1; then
        GMT_BIN=$(which gmt)
        # Get the bin directory, then go up one level to get GMT home
        GMT_BIN_DIR=$(dirname "$GMT_BIN")
        export GMTHOME=$(dirname "$GMT_BIN_DIR")
        export GMT4HOME="$GMTHOME"  # For HC backward compatibility
    elif [ -d "/usr/include/gmt" ] && [ -f "/usr/lib/x86_64-linux-gnu/libgmt.so" ]; then
        # Ubuntu/Debian system GMT6
        export GMTHOME="/usr"
        export GMT4HOME="$GMTHOME"  # For HC backward compatibility
    elif [ -d "/usr/local/include/gmt" ] && [ -f "/usr/local/lib/libgmt.dylib" ]; then
        # macOS Homebrew GMT6
        export GMTHOME="/usr/local"
        export GMT4HOME="$GMTHOME"
    fi
fi
if [ -z "$NETCDFHOME" ]; then
    if [ "$(uname -s)" = "Darwin" ]; then
        if command -v brew >/dev/null 2>&1 && brew ls --versions netcdf-fortran >/dev/null 2>&1; then
            export NETCDFHOME="$(brew --prefix netcdf-fortran)"
            log_info "Detected Homebrew netcdf-fortran: $NETCDFHOME"
            if brew ls --versions netcdf >/dev/null 2>&1; then
                export NETCDF_C_HOME="$(brew --prefix netcdf)"
                log_info "Detected Homebrew netcdf (C library): $NETCDF_C_HOME"
            fi
        else
            prefix=$(detect_system_netcdf) || true
            if [ -n "$prefix" ] && [ -d "$prefix" ]; then
                export NETCDFHOME="$prefix"
                log_info "Detected system NetCDF via nc-config: $NETCDFHOME"
            else
                export NETCDFHOME="$(pwd)/netcdf-c-4.9.3-rc1"
                log_info "Homebrew netcdf-fortran not found; falling back to bundled NetCDF at $NETCDFHOME"
            fi
        fi
    else
        prefix=$(detect_system_netcdf) || true
        if [ -n "$prefix" ] && [ -d "$prefix" ]; then
            export NETCDFHOME="$prefix"
            log_info "Detected system NetCDF via nc-config: $NETCDFHOME"
        else
            export NETCDFHOME=$(pwd)/netcdf-c-4.9.3-rc1
            log_info "Using bundled NetCDF at $NETCDFHOME"
        fi
    fi
else
    log_info "NETCDFHOME preset to $NETCDFHOME (honoring user override)"
fi
export ARCH=$(uname -m)

# Create symlink only if it doesn't exist
if [ ! -e "$NETCDFHOME/lib/libnetcdf.so.7" ]; then
    ln -s $NETCDFHOME/lib/libnetcdf.so.22 $NETCDFHOME/lib/libnetcdf.so.7 2>/dev/null || true
fi

export PATH=$PATH:$NETCDFHOME/lib:$GMT_GSHHG_DATA
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDFHOME/lib
