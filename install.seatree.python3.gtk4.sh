#! /bin/bash

# The shell script is to set up environment for SEATREE.
# Currently supported systems include:
# 1. ubuntu: Ubuntu 22.04
# 2. docker: please see Dockerfile.
# 3. macos: MacOS 14

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
          sudo apt-get install git wget cmake python3 python3-pip python3-numpy
          sudo apt-get install python3-gi python3-gi-cairo gir1.2-gtk-4.0 libgtk-4-dev
          sudo apt-get install build-essential gfortran
          sudo apt-get install x11-apps libx11-dev libxt-dev libxaw7-dev ghostscript libhdf5-dev gawk
          pip3 install --user matplotlib==3.9.2 --break-system-packages # needs newer version of matplotlib to work.
        fi 
    elif [ $MACHINE == "macos" ]; then 
        echo "Installing SEATREE on MacOS 14 ... ..."
        if [ -n "$ENV" ]; then
            brew install gtk4
            pip3 install --break-system-packages matplotlib==3.9.2 pygobject       
        fi   
    fi 

    export SEATREE=$(pwd)
    export GMT4HOME=$(pwd)/gmt-4.5.18
    export GMTHOME=$GMT4HOME
    export NETCDFHOME=$(pwd)/netcdf-c-4.9.3-rc1

    if [ -e "netcdf-c-4.9.3-rc1" ]; then 
        echo "Seems netcdf-c-4.9.3-rc1 has been installed."
    else
        echo "Local netcdf is not available, installing netcdf-c-4.9.3-rc1 ..."
        bash install.netcdf.ubuntu22.sh
    fi 
    if [ -e "gmt-4.5.18" ]; then  
        echo "Seems gmt-4.5.18 has been installed."
    else 
        echo "Local gmt4 is not available, installing gmt4-5.18 ..."
        bash install.gmt4.ubuntu22.sh
    fi 
    yes '' | ./configure.python3.gtk4 

    echo "Installing new ConMan v3.0.0 from CIG ConMan GitHub ..."
    bash install.conman.sh
fi

export SEATREE=$(pwd)
export GMT4HOME=$(pwd)/gmt-4.5.18
export GMTHOME=$GMT4HOME
export GMT_GSHHG_DATA=$GMT4HOME/gshhg-gmt-2.3.7
export NETCDFHOME=$(pwd)/netcdf-c-4.9.3-rc1 
ln -s $NETCDFHOME/lib/libnetcdf.so.22 $NETCDFHOME/lib/libnetcdf.so.7
export PATH=$PATH:$NETCDFHOME/lib:$GMT_GSHHG_DATA
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDFHOME/lib
