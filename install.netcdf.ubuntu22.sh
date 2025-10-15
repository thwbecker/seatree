VERSION="4.9.3-rc1"
SEATREEROOT=$(pwd)
echo "SEATREE ROOT path is" $SEATREEROOT
wget https://github.com/Unidata/netcdf-c/archive/refs/tags/v$VERSION.tar.gz
tar -zxvf v$VERSION.tar.gz
NETCDF_DIR_NAME=netcdf-c-$VERSION
cd $NETCDF_DIR_NAME
mkdir build 
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$SEATREEROOT/$NETCDF_DIR_NAME
make
make install

