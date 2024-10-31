#! /bin/bash

# The shell script is to set up environment for SEATREE.
# Currently supported systems include:
# 1. ubuntu: Ubuntu 22.04
# 2. docker: ubuntu 22.04 with docker.
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
        c)
            MACH=$OPTARG
            CONFIG="True"
            ;;
        h)
            echo "Usage: ./install.seatree.sh [-h] [-m Machine_name] [-c Machine_name] "
            echo "                                                                     "
            echo "Examples:                                                            "
            echo "                                                                     "
            echo "./install.seatree.sh -h                                               "
            echo " -----Display this help message                                      "
            echo "                                                                     "
            echo "./install.seatree.sh -m ubuntu                                           "
            echo " -----Install SEATREE on ubuntu 22.04                            "
            echo "                                                                     "
            echo "./install.seatree.sh -c macos                                        "
            echo " -----Simply set up envs for SEATREE without installation             "
            echo " -----on macos                                                      "
            echo "                                                                     "
            echo "source install.seatree.sh                                             "
            echo " -----Activate ENV VAR EQQUASIROOT and add exes to PATH              "
            echo "                                                                     "
            echo "Currently supported machines include:                                "
            echo " ubuntu/macos/docker                                                    "
            ;;
    esac
done 

if [ -n "$MACH" ]; then 
    export MACHINE=$MACH
    if [ $MACHINE == "ubuntu" ]; then 
        if [ -n "$ENV" ]; then
          sudo apt-get install git vim python3 python3-pip python3-matplotlib gmt python3-numpy
          sudo apt-get install python3-gi python3-gi-cairo gir1.2-gtk-4.0 libgtk-4-dev
          sudo apt-get install build-essential
          sudo apt-get install gfortran
          pip3 install --user matplotlib==3.9.2 --break-system-packages # needs newer version of matplotlib to work.
        fi 
    elif [ $MACHINE == "docker" ]; then 
        echo "Installing SEATREE on Docker Ubuntu 22.04 ... ..."
        if [ -n "$ENV" ]; then
          apt-get install git vim python3 python3-pip python3-matplotlib gmt python3-numpy
          apt-get install python3-gi python3-gi-cairo gir1.2-gtk-4.0 libgtk-4-dev
          apt-get install build-essential
          apt-get install gfortran
          pip3 install --user matplotlib==3.9.2 --break-system-packages # needs newer version of matplotlib to work.
        fi 
    elif [ $MACHINE == "macos" ]; then 
        echo "Installing SEATREE on MacOS 14 ... ..."
        if [ -n "$ENV" ]; then
            brew install gtk4
            pip3 install --break-system-packages matplotlib==3.9.2 pygobject       
        fi   
    fi 

    export GMT4HOME=$(pwd)/gmt-4.5.18
    export GMTHOME=$GMT4HOME
    export NETCDFHOME=$(pwd)/netcdf-c-4.9.3-rc1
    yes '' | ./configure.python3.gtk4 
fi

export GMT4HOME=$(pwd)/gmt-4.5.18
export GMTHOME=$GMT4HOME
export GMT_GSHHG_DATA=$GMT4HOME/gshhg-gmt-2.3.7
export NETCDFHOME=$(pwd)/netcdf-c-4.9.3-rc1 
ln -s $NETCDFHOME/lib/libnetcdf.so.22 $NETCDFHOME/lib/libnetcdf.so.7
export PATH=$PATH:$NETCDFHOME/lib:$GMT_GSHHG_DATA
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDFHOME/lib
