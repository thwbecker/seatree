#! /usr/bin/env python3

import os, sys
import platform

def add_subdirectories_to_path(root_dir):
    for dirpath, dirnames, filenames in os.walk(root_dir):
        sys.path.append(dirpath)

# Replace '/path/to/python' with the path to your 'python' folder
root_directory = os.getcwd()
add_subdirectories_to_path(root_directory)

# Set up SEATREE environment variable
os.environ["SEATREE"] = root_directory

# Set up ARCH environment variable (auto-detect architecture)
os.environ["ARCH"] = platform.machine()

# Set up GMT environment variables for direct execution
gmt_home = os.path.join(root_directory, "gmt-4.5.18")
if os.path.exists(gmt_home):
    os.environ["GMT4HOME"] = gmt_home
    os.environ["GMTHOME"] = gmt_home
    os.environ["GMT_SHAREDIR"] = os.path.join(gmt_home, "share")
    os.environ["GMT_DATADIR"] = os.path.join(gmt_home, "share")
    os.environ["GMT_FONTPATH"] = os.path.join(gmt_home, "share", "pslib")
    os.environ["GMT_GSHHG_DATA"] = os.path.join(gmt_home, "gshhg-gmt-2.3.7")

    # Add GMT to PATH
    gmt_bin = os.path.join(gmt_home, "bin")
    if gmt_bin not in os.environ.get("PATH", ""):
        os.environ["PATH"] = gmt_bin + os.pathsep + os.environ.get("PATH", "")

# Set up NetCDF environment variables if needed
netcdf_home = os.path.join(root_directory, "netcdf-c-4.9.3-rc1")
if os.path.exists(netcdf_home):
    os.environ["NETCDFHOME"] = netcdf_home

    # Add NetCDF lib to LD_LIBRARY_PATH
    netcdf_lib = os.path.join(netcdf_home, "lib")
    if "LD_LIBRARY_PATH" in os.environ:
        os.environ["LD_LIBRARY_PATH"] = netcdf_lib + os.pathsep + os.environ["LD_LIBRARY_PATH"]
    else:
        os.environ["LD_LIBRARY_PATH"] = netcdf_lib

#for path in sys.path:
#    print(path)

os.system('python3 '+root_directory+'/python3/seatree/gui/SEATREE.py')
