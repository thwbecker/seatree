#! /bin/bash

SUDO='sudo'

if [ $1 == 'docker' ]; then
  SUDO=' ' 
fi 

$SUDO apt-get install git vim python3 python3-pip python3-matplotlib gmt python3-numpy
$SUDO apt-get install python3-gi python3-gi-cairo gir1.2-gtk-4.0 libgtk-4-dev
$SUDO apt-get install build-essential
$SUDO apt-get install gfortran
pip3 install --user matplotlib==3.9.2 # needs newer version of matplotlib to work.

