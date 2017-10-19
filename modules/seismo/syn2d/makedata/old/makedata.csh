#!/bin/csh
# example: 
# square region 
set xtot=100.
# equal parameterization increment in both directions:
set dx=1.
### xtot and dx should be consistent with values assigned in makemodel
# number of pixels in parameterization grid ((xtot/dx)^2):
set nfree = 10000 
# root of name of output files with ray geometry info:
set name=rays
# total number of data to be generated:
set ndata=20000
# minimum acceptable epicentral distance:
set deltamin=50.
# increment along ray path:
set rpinc=0.1
# input model:
set model=../makemodel/chkbd.px
# output data file:
set namedata=rays.chkbd.rhs
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
./datamaker<<EOF
$name.xxx
$name.ind
$name.pnt
$model
$namedata
$nfree
$ndata
EOF
