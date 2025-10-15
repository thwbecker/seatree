#!/bin/bash

# GMT 4.5.18 Installation Script for Ubuntu 22.04

# Define installation directory
SEATREEROOT=$(pwd)
GMT_INSTALL_DIR=$SEATREEROOT"/gmt-4.5.18"
GSHHG_DIR=$GMT_INSTALL_DIR"/gshhg-gmt-2.3.7"

rm -rf gmt-4.5.18
if [ -e "gmt-4.5.18-src.tar.bz2" ]; then
    echo "gmt-4.5.18-src.tar.bz2 exists."
else
    echo "Downloading GMT 4.5.18 source code..."
    wget ftp://ftp.soest.hawaii.edu/gmt/gmt-4.5.18-src.tar.bz2
fi

echo "Extracting GMT 4.5.18 source..."
tar -xvf gmt-4.5.18-src.tar.bz2
cd gmt-4.5.18

# Downloading GSHHG dataset
wget http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-gmt-2.3.7.tar.gz
tar -zxvf gshhg-gmt-2.3.7.tar.gz

# Set NETCDFHOME environment variable
export NETCDFHOME=$SEATREEROOT"/netcdf-c-4.9.3-rc1"
echo $NETCDFHOME
# Configure the GMT build
echo "Configuring GMT 4.5.18..."
./configure --prefix=$GMT_INSTALL_DIR --enable-netcdf=$NETCDFHOME --with-gshhg-dir=$GSHHG_DIR

# Compile the GMT source code
echo "Compiling GMT 4.5.18..."
make
# Install GMT 4.5.18
echo "Installing GMT 4.5.18 to $GMT_INSTALL_DIR..."
make install

# Set environment variables
echo "Setting up environment variables..."
echo "export GMT4HOME=$GMT_INSTALL_DIR" >> ~/.bashrc
echo "export PATH=\$GMT_HOME/bin:\$PATH" >> ~/.bashrc
source ~/.bashrc

# Verify the installation
echo "Verifying GMT 4.5.18 installation..."
gmt --version

if [ $? -eq 0 ]; then
    echo "GMT 4.5.18 installed successfully!"
    echo "GMT4HOME is ", $GMT4HOME
else
    echo "There was an issue with the installation."
fi
# 
