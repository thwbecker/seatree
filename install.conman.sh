#! /bin/bash
export SEATREEROOT=$(pwd)
cd modules/mc
mv ConMan ConMan.legacy #rename existing legacy version
# Fetching the latest ConMan from CIG
git clone https://github.com/geodynamics/conman.git ConMan
cd ConMan
cd src
make
cd $SEATREEROOT


