#!/bin/bash

# GMT 4.5.18 Installation Script for Ubuntu 22.04

# Define installation directory
GMT_INSTALL_DIR=$(pwd)"/gmt-4.5.18"

# Update system and install required dependencies
echo "Updating system and installing dependencies..."
#sudo apt update
#sudo apt install -y build-essential gfortran cmake libx11-dev libnetcdf-dev libgdal-dev \
#libpcre3-dev libcurl4-openssl-dev libfftw3-dev liblapack-dev libblas-dev wget

rm -rf gmt-4.5.18
echo "Downloading GMT 4.5.18 source code..."
wget ftp://ftp.soest.hawaii.edu/gmt/gmt-4.5.18-src.tar.bz2
echo "Extracting GMT 4.5.18 source..."
tar -xvf gmt-4.5.18-src.tar.bz2
cd gmt-4.5.18

# Set NETCDF_HOME environment variable
export NETCDF_HOME=/usr

# Configure the GMT build
echo "Configuring GMT 4.5.18..."
./configure --prefix=$GMT_INSTALL_DIR

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
