ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
c
c  Include file for program to generate input decks for use by:
c
c    Standard ConMan
c    Double Diffusive (DD) ConMan
c    SCAM
c    Double Diffusive (DD) SCAM
c    Annulus ConMan
c    ChainMan
c    DefMan
c
c  All cards are not represented here. Those entries listed in data
c    statements may be changed to suit the user's preference although
c    some are not used by ConMan other than as place holders. Those
c    listed as parameters are required by current (Aug. 1992) versions
c    of ConMan to have the given values. All relevant values within the 
c    data statements are overwritten by those read in from existing 
c    input decks (if such a reference is chosen by the user) with the 
c    exception of the Velocity Boundary Condition Flag cards, the 
c    Absolute Velocity cards, and the Absolute Temperature and 
c    Composition cards.
c   
c
ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
c
      parameter (maxmat = 10)
      parameter (maxsuf = 100)
c
c  code versions
      parameter (istd = 1)
      parameter (iddc = 2)
      parameter (iscm = 3)
      parameter (idds = 4)
      parameter (iann = 5)
      parameter (ichn = 6)
      parameter (idef = 7)
c
      dimension  visc(maxmat),    alam(maxmat),   diff(maxmat),
     &          diffb(maxmat),      ra(maxmat),    rab(maxmat),
     &           dmhu(maxmat),   estar(maxmat),   toff(maxmat),
     &          vstar(maxmat),   x2ref(maxmat), sigref(maxmat),
     &         viscut(maxmat),
c
     &            nel(maxsuf),   iside(maxsuf),  fnorm(maxsuf),
     &           ftan(maxsuf),    flux(maxsuf)
c
      logical
     &       lmovie
c  create "geom.movie" file?
c
      character*1
     &           sep
c  file name prefix-suffix separator
c
      character*11
     &            fstat
c  file status (unknown, new)
c
      character*80 
     &            sdeflt,
c  default suffix
     &            rdeflt,
c  default prefix for run file
c
     &            ideflt,
c  default prefix for (new) input deck
     &            gdeflt,
c  default prefix for (new) geometry deck
     &            mdeflt,
c  default prefix for (new) movie geometry deck
c
     &            outfn,
c  default prefix for output file
     &            rinfn,
c  default prefix for input restart file
     &            routfn,
c  default prefix for output restart file
     &            tsfn,
c  default prefix for time series file
     &            fldfn,
c  default prefix for field file
     &            meanfn,
c  default prefix for mean properties file
     &            strfn,
c  default prefix for stress / strain rate field file
     &            cordfn,
c  default prefix for coordinate file (unformatted version only)
c
c  For ChainMan only:
     &            richfn,
c  default prefix for input restart file for chain link location file
     &            rochfn,
c  default prefix for output restart file for chain link location file
     &            chnfn
c  default prefix for chain link location file
c
c  File status may be changed by the user to either 'unknown' or 'new'
c    depending on whether the user is willing to risk mistakenly 
c    overwriting existing files!
      parameter (fstat = 'unknown')
c
c  Default prefix-suffix separator may be changed by the user.
      parameter (sep = '.')
c
      data sdeflt
     &    / 'new'   /
c
      data mdeflt
     &    / 'movie' /
c
      data lmovie
     &    / .true.  /
c
      data rdeflt
     &    / 'run'   /
c
      data ideflt,  gdeflt,   outfn,   rinfn,  routfn,    tsfn,   fldfn,
     &     meanfn,   strfn,  cordfn
     &    /  'in',  'geom',   'out',   'rin',  'rout',  'tser', 'field',
     &     'mean',   'str', 'coord'   /
c
c  For ChainMan only:
      data richfn,  rochfn,   chnfn
     &   /'richn', 'rochn',  'chain'  /
c
c  Code Version
      data icode / istd /
c
c  Title Card
      data ititle / 'A short title goes here' /
c
c  Global Constants Card
c
c  predetermined constants
      parameter (nsd = 2)
      parameter (ndof = 2)
c
c  determined at run time by "gendeck"
      data  nwrap
     &    /     0   /
c
c  unimplemented options
      data  ntseq,   numeg
     &    /     1,       1   /
c
      data   nelx,    nelz,   iflow,   necho,  inrstr,  iorstr,  nstres,
     &     ntimvs,    isky,   lwork,   nnnit,    expo
     &    /    32,      32,       1,       2,       0,       1,       1,
     &          0,       1,       0,       3,    0.25   /
c
c  Lenardic & Kaula [1993] filter Card (for DD ConMan and DD SCAM only)
c
      data ilkflt
     &    /     0   /
c
c  Grid Deformation Parameter Card (for DefMan only)
c
      data   ibnd,  igrdbt,  igngrd,  srfden
     &    /     1,       1,       0,     1.0   /
c
c  Time Sequence Card
c
c  predetermined constants
      parameter (niter = 2)
      parameter (alpha = 0.5)
c
c  unimplemented option
      data epstol
     &    /1.0e-6   /
c
      data  nstep,    delt,   dtfrc
     &    /  1000,       1,       1   /
c
c  Output Step Card
c
c  unimplemented option
      data nsmprt
     &    /    50   /
c
      data nsdprt,  nsvprt,  nstprt
     &    /    50,      50,      50   /
c
c  Velocity Boundary Condition Flag Cards
c
      data   ivxb,    ivzb,    ivxt,    ivzt,    ivxl,    ivzl,    ivxr,
     &       ivzr
     &    /     0,       1,       0,       1,       1,       0,       1,
     &          0   /
c
c  (modified) Initial Temperature/Composition Card
c
c    "xmin", "xmax", "zmin", "zmax" substitute for "xsize" and "zsize"
c
      data  pertt,   pertb,    xmin,    xmax,    zmin,    zmax
     &    /  0.01,     0.0,     0.0,     1.0,     0.0,     1.0   /
c
c  Element Parameter Card
c
c  predetermined constants
      parameter (ntype = 2)
      parameter (nen = 4)
      parameter (nenl = 4)
      parameter (nedof = 2)
      parameter (nitp = 5)
c
c  unimplemented options
      data  implv,   implt
     &    /     0,       0   /
c
      data  numat,  numsuf
     &    /     1,       0   /
c
c  Viscosity Card
      data   (visc(i), i = 1, maxmat)
     &    /   1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,
     &        1.0,     1.0,     1.0   /
c
c  Penalty Card
      data (  alam(i), i = 1, maxmat)
     &    / 1.0e7,   1.0e7,   1.0e7,   1.0e7,   1.0e7,   1.0e7,   1.0e7,
     &      1.0e7,   1.0e7,   1.0e7   /
c
c  Thermal Diffusivity Card
      data   (diff(i), i = 1, maxmat)
     &    /   1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,
     &        1.0,     1.0,     1.0   /
c
c  Compositional Diffusivity Card
      data  (diffb(i), i = 1, maxmat)
     &    /  0.01,    0.01,    0.01,    0.01,    0.01,    0.01,    0.01,
     &       0.01,    0.01,    0.01   /
c
c  Thermal Buoyancy Card
      data     (ra(i), i = 1, maxmat)
     &    / 1.0e5,   1.0e5,   1.0e5,   1.0e5,   1.0e5,   1.0e5,   1.0e5,
     &      1.0e5,   1.0e5,   1.0e5   /
c
c  Compositional Buoyancy Card
      data    (rab(i), i = 1, maxmat)
     &    /-1.0e5,  -1.0e5,  -1.0e5,  -1.0e5,  -1.0e5,  -1.0e5,  -1.0e5,
     &     -1.0e5,  -1.0e5,  -1.0e5   /
c
c  Internal Heating Parameter Card
      data   (dmhu(i), i = 1, maxmat)
     &    /   0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
     &        0.0,     0.0,     0.0   /
c
c  Activation Energy Card
      data  (estar(i), i = 1, maxmat)
     &    /   0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
     &        0.0,     0.0,     0.0   /
c
c  Temperature-Dependent Viscosity Temperature Offset Card
      data   (toff(i), i = 1, maxmat)
     &    /   0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
     &        0.0,     0.0,     0.0   /
c
c  Activation Volume Card
      data  (vstar(i), i = 1, maxmat)
     &    /   0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
     &        0.0,     0.0,     0.0   /
c
c  x2 Reference Card
      data  (x2ref(i), i = 1, maxmat)
     &    /   0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
     &        0.0,     0.0,     0.0   /
c
c  Reference Stress Card (for Stress-Dependent Rheology)
      data (sigref(i), i = 1, maxmat)
     &    /   1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,
     &        1.0,     1.0,     1.0   /
c
c  Viscosity Cut-Off Card (for Temperature-Dependent Rheology)
      data (viscut(i), i = 1, maxmat)
     &    / 1.0e3,   1.0e3,   1.0e3,   1.0e3,   1.0e3,   1.0e3,   1.0e3,
     &      1.0e3,   1.0e3,   1.0e3   /
c
c  Surface Force / Flux Card
      data nel(1),iside(1),fnorm(1), ftan(1), flux(1)
     &    /     1,       1,     1.0,     1.0,     1.0   /
c
c  Absolute Velocity Card
      data vxbotl,  vxbotr,  vzbotl,  vzbotr,  vxtopl,  vxtopr,  vztopl,
     &     vztopr
     &    /   0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
     &        0.0   /
c
c  Absolute Temperature Card
      data   tbot,    ttop
     &    /   1.0,     0.0   /
c
c  Absolute Composition Card
      data   bbot,    btop
     &    /   1.0,     0.0   /
