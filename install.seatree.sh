#! /bin/bash

# The shell script is to set up environment for SEATREE.
# Currently supported systems include:
# 1. ubuntu: Ubuntu 22.04
# 2. docker: please see Dockerfile.
# 3. macos: MacOS 14

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

log_info "=========================================="
log_info "SEATREE Installation Log"
log_info "Started: $(date)"
log_info "Log file: $LOGFILE"
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
            echo "Usage: ./install.seatree.python3.gtk4.sh [-h] [-m Machine_name] [-e Machine_name] "
            echo "                                                                     "
            echo "Examples:                                                            "
            echo "                                                                     "
            echo "./install.seatree.python3.gtk4.sh -h                                 "
            echo " -----Display this help message                                      "
            echo "                                                                     "
            echo "./install.seatree.python3.gtk4.sh -m ubuntu                          "
            echo " -----Install SEATREE on ubuntu 22.04                                "
            echo "                                                                     "
            echo "./install.seatree.python3.gtk4.sh -e macos                           "
            echo " -----Simply set up envs for SEATREE without installation            "
            echo " -----on macos                                                       "
            echo "                                                                     "
            echo "source ./install.seatree.python3.gtk4.sh                             "
            echo " -----Activate ENV VAR and add exes to PATH                          "
            echo "                                                                     "
            echo "Currently supported machines include:                                "
            echo " ubuntu/macos/docker                                                 "
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
          pip3 install --user matplotlib==3.9.2 #--break-system-packages # needs newer version of matplotlib to work.
        fi 
    elif [ $MACHINE == "macos" ]; then 
        echo "Installing SEATREE on MacOS 14 ... ..."
        if [ -n "$ENV" ]; then
            brew install gtk4 ghostscript gawk
            pip3 install --break-system-packages matplotlib==3.9.2 pygobject       
        fi   
    fi 

    export SEATREE=$(pwd)
    export GMT4HOME=$(pwd)/gmt-4.5.18
    export GMTHOME=$GMT4HOME
    export NETCDFHOME=$(pwd)/netcdf-c-4.9.3-rc1
    export ARCH=$(uname -m)

    log_info "[STEP 1/4] $(date +"%Y-%m-%d %H:%M:%S") - Checking NetCDF installation..."
    if [ -e "netcdf-c-4.9.3-rc1" ]; then
        log_info "  -> NetCDF already installed, skipping."
    else
        log_info "  -> Installing netcdf-c-4.9.3-rc1..."
        bash install/install.netcdf.ubuntu22.sh >> "$LOGFILE" 2>&1
        if [ $? -eq 0 ]; then
            log_info "  -> NetCDF installation completed successfully."
        else
            log_info "  -> ERROR: NetCDF installation failed. Check $LOGFILE for details."
        fi
    fi

    log_info "[STEP 2/4] $(date +"%Y-%m-%d %H:%M:%S") - Checking GMT installation..."
    if [ -e "gmt-4.5.18" ]; then
        log_info "  -> GMT already installed, skipping."
    else
        log_info "  -> Installing gmt-4.5.18..."
        bash install/install.gmt4.ubuntu22.sh >> "$LOGFILE" 2>&1
        if [ $? -eq 0 ]; then
            log_info "  -> GMT installation completed successfully."
        else
            log_info "  -> ERROR: GMT installation failed. Check $LOGFILE for details."
        fi
    fi

    log_info "[STEP 3/4] $(date +"%Y-%m-%d %H:%M:%S") - Configuring Python3/GTK4..."
    yes '' | ./install/configure.python3.gtk4 >> "$LOGFILE" 2>&1
    if [ $? -eq 0 ]; then
        log_info "  -> Configuration completed successfully."
    else
        log_info "  -> ERROR: Configuration failed. Check $LOGFILE for details."
    fi

    log_info "[STEP 4/4] $(date +"%Y-%m-%d %H:%M:%S") - Installing ConMan v3.0.0..."
    bash install/install.conman.sh >> "$LOGFILE" 2>&1
    if [ $? -eq 0 ]; then
        log_info "  -> ConMan installation completed successfully."
    else
        log_info "  -> ERROR: ConMan installation failed. Check $LOGFILE for details."
    fi

    log_info ""
    log_info "=========================================="
    log_info "Installation completed: $(date)"
    log_info "=========================================="
    log_info ""
    log_info "Running verification script..."
    python3 verify_seatree.py
fi

export SEATREE=$(pwd)
export GMT4HOME=$(pwd)/gmt-4.5.18
export GMTHOME=$GMT4HOME
export GMT_GSHHG_DATA=$GMT4HOME/gshhg-gmt-2.3.7
export NETCDFHOME=$(pwd)/netcdf-c-4.9.3-rc1
export ARCH=$(uname -m)

# Create symlink only if it doesn't exist
if [ ! -e "$NETCDFHOME/lib/libnetcdf.so.7" ]; then
    ln -s $NETCDFHOME/lib/libnetcdf.so.22 $NETCDFHOME/lib/libnetcdf.so.7
fi

export PATH=$PATH:$NETCDFHOME/lib:$GMT_GSHHG_DATA
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDFHOME/lib
