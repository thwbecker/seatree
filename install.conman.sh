#! /bin/bash
export SEATREEROOT=$(pwd)
cd modules/mc
if [ -e "ConMan.legacy"]; then 
    rm -rf ConMan
else
    mv ConMan ConMan.legacy #rename existing legacy version
fi 
# Fetching the latest ConMan from CIG
git clone https://github.com/geodynamics/conman.git ConMan
cd ConMan
cd src
make
cp ../conman conman.exp # make a copy of exe in the name of conman.exp
cd $SEATREEROOT


