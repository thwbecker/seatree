#!/bin/csh
# example: 
# square region 
set xtot=100.
# equal parameterization increment in both directions:
set dx = 1.
# length (in parameterization pixels) of side of checkerboard elements:
set dcheck = 21
# abs(amplitude of anomalies):
set anoma = 0.5
set scale = ../pm1.cpt
./chkbd<<EOF
$dx
$dx
$xtot
$xtot
$dcheck
$dcheck
$anoma
EOF
./plot2d.csh chkbd.xyz $scale
#ghostview chkbd.xyz.ps &

