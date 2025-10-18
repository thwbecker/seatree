#! /bin/bash 

# The shell script is to set up environments for EQdyna and 
#	install it. It will call the makefile inside src/ and generate 
#	an executable eqdyna and move it to bin/.

# Currently, the machines supported are:
#	ls6:	Lonestar6 at TACC
#	ubuntu: Ubuntu 22.04

# Usage: install-eqdyna.sh [-h] [-m Machine_name] [-c Machine_name]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

while getopts "m:e:c:h" OPTION; do
    case $OPTION in 
        m)
            MACH=$OPTARG
            ;;
        e)
            MACH=$OPTARG
            ENV="True"
            ;;
        c)
            MACH=$OPTARG
            CONFIG="True"
            ;;
        h)
            echo "Usage: ./install-eqdyna.sh [-h] [-m Machine_name] [-c Machine_name] "
            echo "                                                                     "
            echo "Examples:                                                            "
            echo "                                                                     "
            echo "./install-eqdyna.sh -h                                               "
            echo " -----Display this help message                                      "
            echo "                                                                     "
            echo "./install-eqdyna.sh -m ls6                                           "
            echo " -----Install EQdyna on Lonestar6 at TACC                            "
            echo "                                                                     "
            echo "./install-eqdyna.sh -c ubuntu                                        "
            echo " -----Simply set up envs for EQdyna without installation             "
            echo " -----on ubuntu                                                      "
            echo "                                                                     "
            echo "source install-eqdyna.sh                                             "
            echo " -----Activate ENV VAR EQQUASIROOT and add exes to PATH              "
            echo "                                                                     "
            echo "Currently supported machines include:                                "
            echo " ls6/ubuntu/grace                                                    "
            ;;
    esac
done 

if [ -n "$MACH" ]; then 
    export MACHINE=$MACH
    if [ $MACHINE == "ls6" ]; then 
        echo "Installing EQdyna on Lonestar6 at TACC ... ..."
        
        echo "Loading netcdf/4.6.2 module ... ..."
        module load netcdf/4.6.2 
        ml
        
        echo "NETCDF INC and LIB PATH"
        echo $TACC_NETCDF_INC
        echo $TACC_NETCDF_LIB
        
    elif [ $MACHINE == "ubuntu" ]; then 
        echo "Installing EQdyna on Ubuntu 22.04 ... ..."
        export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
        if [ -n "$ENV" ]; then
            # It uses MPICH MPI.
            apt-get install git vim make mpich
            apt-get install libnetcdf-dev libnetcdff-dev 
            apt-get install python3 python3-pip
            pip install numpy netCDF4 matplotlib xarray
            pip install --upgrade numpy
        fi 
    elif [ $MACHINE == "grace" ]; then 
        echo "Installing EQdyna on Grace at TAMU ... ..."
        echo "Loading netcdf module ... ..."
        module load netCDF
        ml
        
        echo "NETCDF INC and LIB PATH"
        echo ${EBROOTNETCDF}/include
        echo ${EBROOTNETCDF}/lib64
    
    elif [ $MACHINE == "macos" ]; then
        echo "Installing EQdyna on MacOS ... ..."
        # Use NETCDFHOME if set, otherwise use Homebrew NetCDF
        if [ -z "$NETCDFHOME" ]; then
            echo "NETCDFHOME not set, using Homebrew NetCDF"
        else
            echo "Using NETCDFHOME: $NETCDFHOME"
        fi
        if [ -n "$ENV" ]; then
            brew install mpich python
            pip3 install --break-system-packages numpy matplotlib xarray
            # Install netCDF4 Python package with SEATREE NetCDF
            echo "Installing Python netCDF4 with SEATREE NetCDF: $NETCDFHOME"
            NETCDF_DIR=$NETCDFHOME USE_NCCONFIG=1 pip3 install --no-binary netCDF4 --break-system-packages netCDF4
        fi
    fi 
    
    if [ -n "$CONFIG" ]; then 
        echo "Simply configure EQdyna without installation ... ..."
    else
        cd "$SCRIPT_DIR/src"
        make
        cd "$SCRIPT_DIR"
        mkdir -p bin
        mv src/eqdyna bin
    fi

    export EQDYNAROOT=$SCRIPT_DIR
    export PATH=$SCRIPT_DIR/bin:$PATH
    export PATH=$SCRIPT_DIR/scripts:$PATH
    
    chmod -R 755 "$SCRIPT_DIR/scripts"
fi

export EQDYNAROOT=$SCRIPT_DIR
export PATH=$SCRIPT_DIR/bin:$PATH
export PATH=$SCRIPT_DIR/scripts:$PATH

echo EQDYNAROOT
echo PATH 
