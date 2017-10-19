#! /usr/bin/env python

import os, signal, sys, traceback
from optparse import OptionParser
from pylab import *


# read model data 
a = load('image2.xyz')          # using pylab function write to numpy array A
n = int(sqrt(a.shape[0]))                # assume square

x=a[:,0].reshape(n,n)
y=a[:,1].reshape(n,n)
z=a[:,2].reshape(n,n)


pcolor(x, y,z,cmap=cm.jet,shading='flat')
#
# add sources and receivers for fun
#
src=load('../makedata/sources.txt')
rec=load('../makedata/receivers.txt')

plot(src[:,0],src[:,1],'r^',rec[:,0],rec[:,1],'bo',)

xlim(1,n);ylim(1,n);
colorbar()
title('input model')
show()
