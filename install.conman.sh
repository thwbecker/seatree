#! /bin/bash
export SEATREEROOT=$(pwd)
cd modules/mc
if [ -e "ConMan.legacy" ]; then 
    rm -rf ConMan.legacy
fi

mv ConMan ConMan.legacy #rename existing legacy version

# Fetching the latest ConMan from CIG
git clone https://github.com/geodynamics/conman.git ConMan
cd ConMan
cd src
make
cp ../conman conman.exp # make a copy of exe in the name of conman.exp
cd $SEATREEROOT


