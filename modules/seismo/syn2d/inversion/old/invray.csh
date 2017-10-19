#!/bin/csh
# square region 
set xtot=100.
# equal parameterization increment in both directions
set dx = 1.
#set damp=100.
set name=../makedata/rays
set namedata=../makedata/rays.chkbd-1.rhs
set ndata=20000
rm lcurve.txt
foreach damp (0. 0.01 0.5 1. 5. 7.5 10. 15. 25. 35. 50. 75. 100. 200. 300. 500. 1000.)
set nameout=maps/rays.chkbd.out-$damp
./invray<<EOF
$name
$namedata
$ndata
$dx
$dx
$xtot
$xtot
$nameout
$damp
EOF
./plot2d.csh $nameout.xyz pm1.cpt
end

