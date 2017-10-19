#!/bin/csh
# example: 
# square region 
set xtot=100.
# equal parameterization increment in both directions
set dx=1.
set nfree = 10000 # depends on values above
set name=rays
set ndata=20000
set deltamin=50.
set rpinc=0.1
set model=../makemodel/chkbd.px
./shootray<<EOF
$dx
$dx
$xtot
$xtot
$name
$ndata
1
$deltamin
$rpinc
EOF
# experiment with different levels of noise
foreach sigma (2 5)
set namedata=rays.chkbd-$sigma.rhs
./noisydatamaker<<EOF
$name.xxx
$name.ind
$name.pnt
$model
$namedata
$nfree
$ndata
$sigma
11
EOF
end
