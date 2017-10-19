#!/bin/csh
# arg 1 = xyz file, arg 2 = colorscale
xyz2grd $1 -G$1.grd -R0/100/0/100 -I1/1 -N0 -V 
grdimage $1.grd -R -JX16in  -Y3in -C$2 -K -B10/10WESN -P >! $1.ps
psscale -C$2 -O -D8in/-1in/12in/1inh -V >>$1.ps
