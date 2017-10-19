Generate synthetic data from an input model, and invert them to see how well
the input model is reproduced. 

by Lapo Boschi, Sep 2007;modifications by Thorsten Becker, 2008

(1) First of all, generate an input model.

directory makemodel

run program chkbd (simply compiled by f77 -o chkbd chkbd.f). use script
chkbd.csh to start. This script is particularly simple, the model is defined
on a square (rather than a rectangle, which would also be possible) and
checkerboard elements are also squares. The specs of the input model to be
generated are set at the beginning of chkbd.csh. The output model is generated
in three forms:
chkbd.px  - model value vs. pixel index
chkbd.xyz  - model value vs. x,y location of pixel center
chkbd.xyz.grd  - *.xyz file converted to gmt .grd format
and as a postscript file that can be checked:
chkbd.xyz.ps

Alternatively, use convert_pgm to convert a gray scale (PGM) image
into a model file


(2) Then, generate a synthetic database associated with the model

directory makedata

(2.1)  generate a set of randomly located "sources" and "receivers", then for
each source-receiver couple, "trace a ray" through the model. No raytracing is
actually conducted, but energy  is assumed to propagate along straight lines
(i.e. ray paths are straight lines). For each source-station couple, one still
needs to compute the phase anomaly accumulated along the ray path.

the program shootray.f randomly generates source locations (then saved in file
sources.txt) and, for each random source location, generated a random receiver
location (receivers.txt), such that source-station distance is above the
prescribed threshold (set in the scripts makedata.csh and
makenoisydata.csh). For each source-station couple, the length of the
corresponding ray in each parameterization pixel is also calculated by
shootray.f, and the results are stored in a tomographic "A" matrix, in the
usual sparse format (.xxx, .ind, .pnt files). 

(2.2) the programs datamaker.f and noisedatamaeker.f then dot-multiply the "A"
matrix with a chosen input model (generated as described at step 1 above), to
generate a synthetic dataset, file *.rhs. The only difference between the two
programs, is that one adds random noise to the data. 

the sequence shootray/datamaker or noisedatamaker can be executed using the
corresponding scripts makedata.csh and makenoisydata.csh. Parameters including
input model, total number of data, level of artificial noise, are explained
and can be set in those scripts.

(3) ...We are now ready to invert the synthetic data.

directory inversion

(3.1) the program invray.f finds the least-squares solution associated with the
synthetic data of step 2 above. It should be run using the script invray.csh,
where the locations of "A" matrix files:
set name=../makedata/rays
...and of the synthetic data vector:
set namedata=../makedata/rays.chkbd-1.rhs
are specified. 
Again, parameters xtot and dx set here should have the same values as at steps
1 and 2--they define the parameterization. 
invray.csh runs invray.f repeatedly, experimenting with a range of values for
the damping parameter: 
foreach damp (0. 0.01 0.5 1. 5. 7.5 10. 15. 25. 35. 50. 75. 100. 200. 300. 500. 1000.)
you will therefore end up with a corresponding range of differently damped
least-squares solutions, all placed in subdirectory maps/ in the current
version of invray.csh:
-rw-r--r--  1 boschil users 280000 2007-09-24 18:27 maps/rays.chkbd.out-1.
-rw-r--r--  1 boschil users 280000 2007-09-24 18:27 maps/rays.chkbd.out-5.
-rw-r--r--  1 boschil users 280000 2007-09-24 18:28 maps/rays.chkbd.out-10.
-rw-r--r--  1 boschil users 280000 2007-09-24 18:28 maps/rays.chkbd.out-15.
-rw-r--r--  1 boschil users 280000 2007-09-24 18:28 maps/rays.chkbd.out-25.
-rw-r--r--  1 boschil users 280000 2007-09-24 18:28 maps/rays.chkbd.out-35.
etc...
(plus the corresponding xyz, grd and ps files, and log files keeping track of
number of iterations required for convergence, and corresponding achieved
variance reduction)

Use invray.csh to run invray.f

A by-product of invray.f is the "hitcount" associated with the chosen
parameterization, i.e. the number of rays crossing each block. This is output
as a file called hitcount.xyz (or hitcount.dat, with the usual naming
convention). You can easily plot the xyz file, using plot2d.csh, with an
appropriate GMT colorpalette file (e.g., count550.cpt, which is provided):
./plot2d.csh hitcount.xyz count550.cpt

plot2d.csh is the same script used above to plot solution models,
automatically executed for each solution by invray.csh.

(3.2) 
to evaluate the trade-off between damping (or complexity of the solution) and
datafit (reduction of variance), we can conduct a trade-off or L-curve
analysis. A plot of misfit (1.-variance reduction; this is what is typically
used in the L-curve analysis) vs. model complexity (here, I use model RMS(*),
but this is an arbitrary choice, one could use other functions to quantify the
complexity of a solution) is provided by invray.f: lcurve.txt. 

In the sub-directory curves/, a few simple Fortran programs, executed in the
right sequence by the script curv.csh, compute the "curvature" of the
l-curve. For example (from directory curves):
./curv.csh ../lcurve.txt
now the file lcurve.txt.der.02.wrt.01.der.03.wrt.01.crv
has model complexity on the first column, misfit on the second column, first
derivative of misfit with respect to complexity on the third column, second
derivative on the fourth column, and curvature (a function of all these
quantities) on the fifth (last) column. The point of maximum curvature
coincides with the corner of the L-curve, and some authors argue that this is
the best compromise between exceedingly complex solutions that fit the data
very well (but also fit the noise), and simpler, smoother solutions that do
not explain the data so well (the smoother the solution, the lower the fit).  


Minor modifications by Thorsten Becker, Sep 30, 2007


global_para.dat file determines all settings, pm1.cpt colorbar all
plot styles

csh scripts were converted to bash


run "make"

to automatically generate all steps above
