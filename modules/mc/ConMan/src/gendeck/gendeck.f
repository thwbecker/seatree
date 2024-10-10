ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
c
c  Program to generate input decks for use by:
c
c    Standard ConMan
c    Double Diffusive (DD) ConMan
c    SCAM
c    Double Diffusive (DD) SCAM
c    Annulus ConMan
c    ChainMan
c    DefMan
c
c  Program documentation assumes that the reader is familiar with 
c    ConMan and with the ConMan user's guide.
c
c  This program is not yet smart enough to allow the user to create the 
c    full range of input decks nor will it always prohibit the user from
c    responding inappropriately. Rather, it will facilitate the making of
c    fairly simple input decks which can then be modified by the user
c    using a text editor.
c
c  The user may select (via a carriage return) defaults displayed in
c    [square brackets]. These default values can come from three 
c    sources: an existing input deck, an include file 
c    ('deck.defaults.h'), or from "gendeck" itself.
c
c    If one chooses not to use an existing input deck as a reference,
c      all relevant default values given in 'deck.defaults.h' will be
c      offered to the user. If one chooses to use an existing input deck
c      as a reference, all of these default values except for those
c      associated with the Velocity Boundary Condition Flag cards, the
c      Absolute Velocity cards or the Absolute Temperature and 
c      Composition cards will be replaced by those found within the 
c      reference input deck.
c      Note: text found in reference input decks after the Title Card 
c      may result in a program crash. A subsequent version of "gendeck"
c      may allow and even produce such text.
c
c    The defaults corresponding to the following options are hardwired
c      in "gendeck" and do not depend at all on values present in either
c      'deck.defaults.h' or a selected reference input deck:
c
c        i1. Velocity Boundary Condition Flag Cards
c        i2. Temperature/Composition Boundary Condition Flag Cards
c        i3. Surface Force / Flux Cards
c
c        g1. Coordinate Group Cards
c        g2. Velocity Boundary Condition Group Cards
c        g3. Temperature/Composition Boundary Condition Group Cards
c        g4. Element Connectivity Generation Group Cards
c
c    Filename conventions: GenDeck will create the "in", "geom", and 
c      "run" files with one user-determined suffix. In addition, GenDeck
c      will create a "geom.movie" file for use by "MTGP".  It will also 
c      generate a "run" file with all entries having this same suffix. 
c      Prefixes for all ConMan input and output files may be changed in
c      'deck.defaults.h'.
c
c  1.0  => Scott King (July 1992)
c  2.0  => Steven S. Shapiro (August 1992)
c  2.01 => Steven S. Shapiro (February 1993)
c  2.02 => Steven S. Shapiro (April 1993)
c  2.03 => Steven S. Shapiro (June 1993)
c  2.04 => Steven S. Shapiro (October 1993)
c  2.05 => Steven S. Shapiro (November 1993)
c
ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      program gendeck
c
c
c   ifort didn't think that xmin,xmax,zmin,zmax and such were real
c
      implicit real(a-h, o-z), integer(i-n)

      real*8 pi
c

      character*56 ititle
c  title excluding date and time

      include 'deck.defaults.h'
c
      character*80 ifile, gfile, mfile, rfile, irefnm, suffix, fname 
c     &           , grefnm
      character*1 dfault
c
      logical yes, iref, gref, done
c
      integer corner, getint, gunit, munit, runit, grefu, velbcf
c
      logical lexist
      logical addcomment

      parameter (istdi = 5)
      parameter (istdo = 6)
      parameter (iunit = 10)
      parameter (gunit = 11)
      parameter (munit = 12)
      parameter (runit = 13)
      parameter (irefu = 14)
      parameter (grefu = 15)
c
      pi = acos(-1.d0)

      addcomment = .true.
c
c... banner
c
      write (istdo,"(///,26x,'Welcome to GenDeck',///,
     &               17x,'The interactive input deck generator',/,
     &               30x,'for ConMan',///,
     &               21x,'=====>  Version  2.05  <=====',//)")
c
      write (istdo,"(///,5x,'A <return> at any prompt selects the',
     &         /,1x,'default value given in [square brackets]. Some',
     &         /,1x,'default values can be read from an existing',
     &         /,1x,'input deck - others are preset either in',
     &         /,1x,'deck.defaults.h or within the code itself.',
     &         /,1x,'See comments at the top of the gendeck source', 
     &         /,1x,'code. For more creative input decks one needs to ',
     &         /,1x,'edit the input decks directly. Some options are ',
     &         /,1x,'only available with particular versions of the ',
     &         /,1x,'family of ConMan codes. Consult the ConMan README',
     &         /,1x,'file for up-to-date information.',
     &        //,1x,'Warning: GenDeck does not (yet!) check for ',
     &         /,1x,'inconsistant or inappropriate entries.',///)")
c
      write (istdo, *)
      write (istdo, *) istd, '=> Standard ConMan'
      write (istdo, *) iddc, '=> Double-Diffusive ConMan'
      write (istdo, *) iscm, '=> SCAM'
      write (istdo, *) idds, '=> Double-Diffusive SCAM'
      write (istdo, *) iann, '=> Annulus ConMan'
      write (istdo, *) ichn, '=> ChainMan'
      write (istdo, *) idef, '=> DefMan'
      icode = getint ('Enter ConMan version:', icode)
c
c... option to open existing "in" and "geom" to use as reference files
c
      write (istdo, *)
      write (istdo, *) 'Reference files must correspond to the version o
     &f ConMan selected above.'
      write (istdo, *)
5     write (istdo, "('Enter name of reference IN file if any:   ', $)")
      read (istdi, '(a)') irefnm
      if (irefnm(1:1) .eq. ' ') then
        iref = .false.
      else
        if (lexist (irefnm)) then
          iref = .true.
          open (unit=irefu, file=irefnm, status='old')
        else
          write (istdo,"(/,'>>>>>>>> Error opening [',a,'] ')")
     &           irefnm (1:kblnk(irefnm))
          write (istdo, *)
          go to 5
        end if
      end if
c
c6     write (istdo, "('Enter name of reference GEOM file if any: ', $)")
c      read (istdi, '(a)') grefnm
c      if (grefnm(1:1) .eq. ' ') then
        gref = .false.
c      else
c        if (lexist (grefnm)) then
c          gref = .true.
c          open (unit=grefu, file=grefnm, status='old')
c        else
c          write (istdo,"(/,'>>>>>>>> Error opening [',a,'] ')")
c     &           grefnm (1:kblnk(grefnm))
c          write (istdo, *)
c          go to 6
c        end if
c      end if
c
c... open "in", "geom", "geom.movie", and "run" files
c
      write (istdo, *)
      write (istdo, *) 'File names are limited to 80 characters in ConMa
     &n.'
      write (istdo, *) 'Warning: GenDeck does not check filename length!
     &'
      write (istdo, *)
c
 10   write (istdo,"('Enter suffix for input filenames [',a,'] ',$)") 
     &       sdeflt (1:kblnk(sdeflt))
      read (istdi,"(a80)") suffix
c
      if (suffix(1:1) .eq. ' ') then 
        suffix = sdeflt
      else
      end if
c
      ifile = ideflt (1:kblnk(ideflt)) // sep // suffix
      gfile = gdeflt (1:kblnk(gdeflt)) // sep // suffix
      mfile = gdeflt (1:kblnk(gdeflt)) // sep // 
     &        mdeflt (1:kblnk(mdeflt)) // sep // suffix
      rfile = rdeflt (1:kblnk(rdeflt)) // sep // suffix
c
      open (unit=iunit, file=ifile, err=11)
      go to 12
11    write (istdo,"(/,'>>>>>>>> Error opening [',a,'] ')")
     &       ifile (1:kblnk(ifile))
      write (istdo, *)
      go to 10
12    continue
c
      open (unit=gunit, file=gfile, err=21)
      go to 22
21    write (istdo,"(/,'>>>>>>>> Error opening [',a,'] ')") 
     &       gfile (1:kblnk(gfile))
      write (istdo, *)
      go to 10
22    continue
c
      if (lmovie) then
        open (unit=munit, file=mfile, err=23)
        go to 24
23      write (istdo,"(/,'>>>>>>>> Error opening [',a,'] ')") 
     &         mfile (1:kblnk(mfile))
        write (istdo, *)
        go to 10
24      continue
      else
      end if
c
      open (unit=runit, file=rfile, err=31)
      go to 32
31    write (istdo,"(/,'>>>>>>>> Error opening [',a,'] ')") 
     &       rfile (1:kblnk(rfile))
      write (istdo, *)
      go to 10
32    continue
c
c... Title Card
c
      if (iref) then
        read (irefu, "(a56)") ititle
      else
c
c... offer default given in 'deck.defaults.h'
c
      end if
      write (istdo, *)
      write (*,"('Enter a title <= 56 characters [',a56,'] ',$)") ititle
      read (istdi, "(a56)") fname
      if (fname(1:1) .ne. ' ') ititle = fname
c
      write (iunit, "(a56)") ititle
c
c... Global Constants Card
c
      if (iref) then
c
c... ignore "nsd", "ndof" - given as parameters in 'deck.defaults.h'
c
        read (irefu, *) numnp, nsdx, ndofx, nelx, nelz, iflow, necho,
     &                  inrstr, iorstr, nstres, nodebn, ntimvs, ntseq,
     &                  numeg, isky, lwork, nwrap, nnnit, expo

      else
c
c... offer defaults given in 'deck.defaults.h'
c
      end if
c
      write (istdo, *)
      write (istdo, *) 'GenDeck assumes a 2-d Cartesian grid with'
      write (istdo, *) '           x1 - horizontal = x'
      write (istdo, *) '           x2 -  vertical  = z'
c
      nelx = getint ('Enter number of x1 elements:', nelx)
      nelz = getint ('Enter number of x2 elements:', nelz)
c
      numel = nelx * nelz 
      numnp = (nelx+1)*(nelz+1)
c
      write (istdo, *)
      write (istdo, *) '0 => terse without initial field data'
      write (istdo, *) '1 => verbose with initial field data'
      write (istdo, *) '2 => terse with initial field data'
      necho = getint ('Enter output option:', necho)
c
      if (iflow .eq. 1) then
        dfault = 'y'
      else
        dfault = 'n'
      end if
c
      write (istdo, *)
      if (yes ('Execute code', dfault)) then
        iflow = 1
      else
        iflow = 0
      end if
c 
      write (istdo, *)
      write (istdo, *) '0 => conductive'
      write (istdo, *) '1 => read from restart file'
      if ((icode .eq. istd) .or. (icode .eq. ichn) .or. 
     &    (icode .eq. idef)) then
        write (istdo, *) '2 => apply boundary layer theory'
      else
      end if
      write (istdo, *) '? => user specific'
      inrstr = getint ('Enter initial buoyancy field option:', inrstr)

      if (iorstr .eq. 1) then
        dfault = 'y'
      else
        dfault = 'n'
      end if

c
      write (istdo, *)
      if (yes ('Write restart file', dfault)) then
        iorstr = 1
      else
        iorstr = 0
      end if
c 
c      write (istdo, *)
c      write (istdo, *) '0 => no stress/strain rate, viscosity, nor effec
c     &tive viscosity output'
c      if ((icode .ne. iscm) .and. (icode .ne. idds)) then
c        write (istdo, *) '1 => stress, viscosity, effective viscosity ou
c     &tput'
c        write (istdo, *) '2 => strain rate, viscosity, effective viscosi
c     &ty output'
c        write (istdo, *) '3 => effective viscosity only'
c      else
c        write (istdo, *) '1 => stress, viscosity output'
c      end if
c      nstres = getint ('Enter stress option:', nstres)
c
      nodebn = 2 * (nelx + 1)
c
      if (ntimvs .eq. 1) then
        dfault = 'y'
      else
        dfault = 'n'
      end if
c
      write (istdo, *)
      if (yes ('Factor stiffness matrix more than once', dfault)) then
        ntimvs = 1
      else
        ntimvs = 0
      end if
c
      if (nwrap .eq. 0) then
        dfault = 'n'
      else
        dfault = 'y'
      end if
c
      if ((icode .ne. iscm) .and. (icode .ne. idds)) then
        write (istdo, *)
        if (yes ('Use wrap around boundary conditions', dfault)) then
          nwrap = nelz
        else
          nwrap = 0
        end if
      else
        nwrap = 0
      end if
c
      if (nwrap .eq. 0) then
        write (istdo, *)
        write (istdo, *) '0 => banded solver'
        write (istdo, *) '1 => skyline solver'
c        write (istdo, *) '2 => dmf solver'
        isky = getint ('Enter solver option:', isky)
c
c        if (isky .eq. 2) then
c          write (istdo, *)
c          lwork = getint ('Enter size of dmf workspace (lwork):', lwork)
c        else
c          lwork = 0
c        end if
c
      else
c
c... can only use skyline solver with wrap around boundary conditions
c
        isky  = 1
        lwork = 0
c
      end if
c
      if ((ntimvs .eq. 1) .and. 
     &    (icode .ne. iscm) .and. (icode .ne. idds) 
     &                      .and. (icode .ne. iann)) then
        write (istdo, *)
c        nnnit = getint ('Enter number of non-Newtonian iterations (>= 1)
c     &:                        ', nnnit)
c        if (nnnit .gt. 1) then
c          expo = gtreal ('Enter weighting of strain-rate- vs. stress-def
c     &fined effective viscosity:', expo)
c        else
c
c... Newtonian rheology
cc
c          expo = 0.0
c
c        end if
           nnnit = 1
          expo  = 0.0
       else
c
c... cannot apply sress-dependence without factoring stiffness matrix 
c      at every time step
c
          nnnit = 1
          expo  = 0.0
c
      end if
c precision
      mprec=2
      if(addcomment)write (iunit,*) 'geometry parameters'
      write (iunit, 9000) numnp, nsd, ndof, nelx, nelz, mprec, iflow, 
     &       necho,
     &       inrstr, iorstr, nodebn, ntimvs, ntseq, numeg, isky,nwrap
c     &       lwork, nwrap, nnnit, expo 
c
      if ((icode .eq. iddc) .or. (icode .eq. idds)) then
        if (iref) then
          read (irefu, *) ilkflt 
        else
        end if
c
        if (ilkflt .eq. 1) then
          dfault = 'y'
        else
          dfault = 'n'
        end if
c
        write (istdo, *)
        if (yes ('Use Lenardic & Kaula (1993) filter', dfault)) then
          ilkflt = 1
        else
          ilkflt = 0
        end if
c
        write (iunit, "(i3)") ilkflt
c
      else if (icode .eq. ichn) then
c
        if (iref) then
          read (irefu, *) nmark 
        else
          nmark = 4 * nelx
        end if
c
        write (istdo, *)
        nmark = getint ('Enter number of marker particles in marker chai
     &n:', nmark)
c
        write (iunit, 9000) nmark
c
      else if (icode .eq. idef) then
c
        if (iref) then
          read (irefu, *) ibnd, igrdbt, igngrd, srfden 
        else


           
c
c... offer defaults given in 'deck.defaults.h'
c
        end if
c

        write (istdo, *)
        igrdbt = getint ('Enter number of nodes (from bottom) to the def
     &orming region:  ', igrdbt)
c
        if (igrdbt .ne. 0) then
c
          ibnd = getint ('Enter number of nodes (from bottom) to the mat
     &erial interface:', ibnd)
c
          if (igngrd .eq. 1) then
            dfault = 'y'
          else
            dfault = 'n'
          end if
          if (yes ('Use routine for initially deformed grid', dfault)) 
     &    then
            igngrd = 1
          else
            igngrd = 0
          end if
c
          srfden = gtreal ('Enter topography constant:', srfden)
c
        else
          ibnd   = 0
          igngrd = 0
          srfden = 0.0
        end if
c
        write (iunit, 9025) igrdbt, ibnd, igngrd, srfden
c
      else
      end if             

c
c... Time Sequence Card
c
      if(addcomment)write (iunit,*) 'time step information'
      if (iref) then
c
c... ignore "niter", "alpha" - given as parameters in 'deck.defaults.h'
c
        read (irefu, *) nstep, niterx, alphax, delt, epstol, dtfrc
c
      else
        if (((icode .eq. iscm) .or. (icode .eq. idds)) .and. 
     &      (dtfrc .eq. 1.0)) then
          dtfrc = 2.0
        else
c
c... offer defaults given in 'deck.defaults.h'
c
        end if
      end if
c
      write (istdo, *)
      nstep = getint ('Enter number of time steps:                 ', 
     &                 nstep)
c
      delt = gtreal ('Enter non-dimensional time to print results:',
     &                delt)
c
      dtfrc = gtreal ('Enter fraction of dt:                       ',
     &                 dtfrc)

      write (iunit, "(i6,i3,f10.6,3(1pe12.5))") nstep, niter, alpha,
     &                                          delt, epstol, dtfrc
c
c... Output Step Card
c
      if (iref) then
        read (irefu, *) nsdprt, nsvprt, nstprt, nsmprt
      else
c
c... offer defaults given in 'deck.defaults.h'
c
      end if
c
      write (istdo, *)
      nsdprt = getint ('Enter steps between restart dumps:      ', 
     &                  nsdprt)
      nstprt = getint ('Enter steps between velocity and temp outputs:',
     &                  nstprt)
      nsmprt = getint ('Enter steps between stress field outputs:     ',   
     &                  nsmprt)
c
      nsvprt = nstprt
      if(addcomment)write (iunit,*) 'output information'
      write (iunit, "(4i6)") nsdprt, nsvprt, nstprt, nsmprt
c
c... Velocity Boundary Condition Flag Cards
c
      write (istdo, *)
      write (istdo, *) '************************************************
     &****************'
      write (istdo, *) 'Velocity Boundary Condition Flag Cards'
      write (istdo, *) 'Only some common options supported.'
      write (istdo, *) 'Not (yet) based on existing input deck or includ
     &e file settings!'
      write (istdo, *) '************************************************
     &****************'
c
      if(addcomment)write (iunit,*) 'velocity boundary conditions'
      if (iref) then
c
c... skip velocity boundary condition flag cards
c
        call skpsec (irefu, 5)
      else
      end if
c
      ifre = 0
      ifix = 1
c
      incx = nelz + 1
      incz = 1
c
c... node numbers cooresponding to the four corners: bottom left,
c      top left, bottom right, top right
c
      ncorbl = 1
      ncortl = nelz + 1
      ncorbr = (incx * nelx) + 1
      ncortr = (nelx + 1) * (nelz + 1)
c
c... bottom edge
c
      write (istdo, *)
      write (istdo, *) '***     B o t t o m    E d g e     ***'
      write (istdo, *) 'Supported options: Unconstrained, Pinned.'
      ivxb = velbcf ('V1=Vx', ivxb, ifre, ifix)
      ivzb = velbcf ('V2=Vz', ivzb, ifre, ifix)
      write (iunit, 9050) ncorbl, ncorbr, incx, ivxb, ivzb
c
c... top edge
c
      write (istdo, *)
      write (istdo, *) '***       T o p      E d g e       ***'
      write (istdo, *) 'Supported options: Unconstrained, Pinned.'
      ivxt = velbcf ('V1=Vx', ivxt, ifre, ifix)
      ivzt = velbcf ('V2=Vz', ivzt, ifre, ifix)
      write (iunit, 9050) ncortl, ncortr, incx, ivxt, ivzt
c
c... left edge
c
      write (istdo, *)
      write (istdo, *) '***       L e f t    E d g e       ***'
      if (nwrap .eq. 0) then
        write (istdo, *) 'Supported options: Unconstrained, Pinned.'
        ivxl = velbcf ('V1=Vx', ivxl, ifre, ifix)
        ivzl = velbcf ('V2=Vz', ivzl, ifre, ifix)
        write (iunit, 9050) ncorbl, ncortl, 1, ivxl, ivzl
      else
        write (istdo, *) 'Wrap around velocity boundary conditions.'
      end if
c
c... right edge
c
      write (istdo, *)
      write (istdo, *) '***      R i g h t    E d g e      ***'
      if (nwrap .eq. 0) then
        write (istdo, *) 'Supported options: Unconstrained, Pinned.'
        ivxr = velbcf ('V1=Vx', ivxr, ifre, ifix)
        ivzr = velbcf ('V2=Vz', ivzr, ifre, ifix)
        write (iunit, 9050) ncorbr, ncortr, 1, ivxr, ivzr
      else
        write (istdo, *) 'Wrap around velocity boundary conditions.'
      end if
c
c... handle corners
c
      if (nwrap .eq. 0) then
        ivxbl = corner (ivxb, ivxl, ifre, ifix)
        ivzbl = corner (ivzb, ivzl, ifre, ifix)
        write (iunit, 9050) ncorbl, ncorbl, 1, ivxbl, ivzbl
c
        ivxbr = corner (ivxb, ivxr, ifre, ifix)
        ivzbr = corner (ivzb, ivzr, ifre, ifix)
        write (iunit, 9050) ncorbr, ncorbr, 1, ivxbr, ivzbr
c
        ivxtl = corner (ivxt, ivxl, ifre, ifix)
        ivztl = corner (ivzt, ivzl, ifre, ifix)
        write (iunit, 9050) ncortl, ncortl, 1, ivxtl, ivztl
c
        ivxtr = corner (ivxt, ivxr, ifre, ifix)
        ivztr = corner (ivzt, ivzr, ifre, ifix)
        write (iunit, 9050) ncortr, ncortr, 1, ivxtr, ivztr
      else
c
c... Arbitrarily fix the bottom left and right corners
c
        write (iunit, 9050) ncorbl, ncorbl, 1, ifix, ifix
        write (iunit, 9050) ncorbr, ncorbr, 1, ifix, ifix
      end if
c
      write (iunit, 9050) 0, 0, 0, 0, 0


c
c... Temperature/Composition Boundary Condition Flag Cards
c

      if(addcomment)write (iunit,*) 'temperature boundary conditions'
      write (istdo, *)
      write (istdo, *) '************************************************
     &**********************'
      write (istdo, *) 'Temperature/Composition Boundary Condition Flag
     &Cards'
      write (istdo, *) 'Only fixed temperature/composition along the top
     & and bottom supported.'
      write (istdo, *) 'Not (yet) based on existing input deck or includ
     &e file settings!'
      write (istdo, *) 'Continue to next card.'
      write (istdo, *) '************************************************
     &**********************'
c
      if ((icode .eq. iddc) .or. (icode .eq. idds)) then
c
c... need to include compositional as well as a thermal field
c
        nfield = 2
      else
        nfield = 1
      end if
c
      do 100 loop = 1, nfield
c
        if (iref) then
c
c... skip temperature/composition boundary condition flag cards
c
          call skpsec (irefu, 4)
        else
        end if
c
c... make bottom and top fixed temperature/composition boundary 
c      conditions
c
c... bottom edge
c
        write (iunit, 9050) ncorbl, ncorbr, incx, ifix
c
c... top edge
c
        write (iunit, 9050) ncortl, ncortr, incx, ifix
c
        write (iunit, 9050) 0, 0, 0, 0
c
100   continue
c
c... Nusselt Number Boundary Condition Flag Cards
c
      if (iref) then
c
c... skip Nusselt Number boundary condition flag cards
c
        call skpsec (irefu, 3)
        call skpsec (irefu, 3)
      else
      end if
c
      write (iunit, 9050) ncorbl, ncorbr, incx
c
      write (iunit, 9050) ncortl, ncortr, incx
c
      write (iunit, 9050) 0, 0, 0
c
      nstart = ncorbl + incz
      nstop  = ncorbr + incz
      write (iunit, 9050) nstart, nstop, incx
c
      nstart = ncortl - incz
      nstop  = ncortr - incz
      write (iunit, 9050) nstart, nstop, incx
c
      write (iunit, 9050) 0, 0, 0
c
c... Initial Temperature/Composition Card
c
      if (iref) then
        if (icode .eq. iddc) then
          read (irefu, *) pertt, pertb, xsize, zsize
        else if (icode .eq. iscm) then
          read (irefu, *) pertt, xsize, zmin, zmax
          zsize = zmax - zmin
        else if (icode .eq. idds) then
          read (irefu, *) pertt, pertb, xsize, zmin, zmax
          zsize = zmax - zmin
        else
          read (irefu, *) pertt, xsize, zsize
        end if
c
      else
c
c... offer defaults based on those given in 'deck.defaults.h'
c
        xsize = xmax - xmin
        zsize = zmax - zmin
      end if
c
      if(addcomment)write (iunit,*) 'initial temp cond information'
    
      write (istdo, *)
      if (inrstr .ne. 1) then
c
c... not using restart file
c
        pertt = gtreal ('Enter temperature perturbation: ', pertt)
      else
      end if
c
      if (((icode .eq. iddc) .or. (icode .eq. idds)) .and. 
     &    (inrstr .ne. 1)) then
c
c... composition code and not using restart file
c
        pertb = gtreal ('Enter composition perturbation: ', pertb)
      else
      end if
c
      if ((icode .eq. iscm) .or. (icode .eq. idds)) then
        xmin  = gtreal ('Enter multiplier of pi for xmin:', xmin)
        xmax  = gtreal ('Enter multiplier of pi for xmax:', xmax)
        xmin  = pi * xmin
        xmax  = pi * xmax
      else
        xmin  = gtreal ('Enter x1 minimum value:         ', xmin)
        xmax  = gtreal ('Enter x1 maximum value:         ', xsize+xmin)
      end if
      zmin  = gtreal ('Enter x2 minimum value:         ', zmin)
      zmax  = gtreal ('Enter x2 maximum value:         ', zsize+zmin)
c
      if (((icode .eq. istd) .or. (icode .eq. ichn) .or. 
     &    (icode .eq. idef)) .and. (inrstr .eq. 2)) then
        rabndy = gtreal ('Enter thermal Rayleigh number to define bounda
     &ry layer thickness:', ra(1))
      else
      end if
c
      if (icode .eq. iddc) then
        write (iunit, "(4f12.6)") pertt, pertb, xmax - xmin, zmax - zmin
      else if (icode .eq. iscm) then
        write (iunit, "(4f12.6)") pertt, xmax - xmin, zmin, zmax
      else if (icode .eq. idds) then
        write (iunit, "(5f12.6)") pertt, pertb, xmax - xmin, zmin, zmax 
      else if (((icode .eq. istd) .or. (icode .eq. ichn) .or. 
     &         (icode .eq. idef)) .and. (inrstr .eq. 2)) then
        write (iunit, "(3f12.6, 1pe12.5)") pertt, xmax - xmin,
     &                                     zmax - zmin, rabndy
      else
        write (iunit, "(4f12.6)") pertt, xmax - xmin, zmax - zmin
      end if
c
c... Element Parameter Card
c
      if(addcomment)write (iunit,*) 'element parameters'
      if (iref) then
c
c... ignore "ntype", "nen", "nenl", "nedof", and "nitp" - given as 
c      parameters in 'deck.defaults.h'
c
c... ignore "numel" - determined from previously entered information
c 
        read (irefu, *) ntypex, numelx, nenx, nenlx, numat, nedofx,
     &                  numsuf, nitpx, implv, implt
      else
c
c... offer defaults given in 'deck.defaults.h'
c
      end if
c
c... material properties
c
      if (iref) then
        if (numat .gt. maxmat) then
          write (istdo, *)
          write (istdo, *) '********************************************
     &********************'
          write (isdto, *) 'Warning: More than the (arbitrary) limit of
     &maxmat = ', maxmat
          write (isdto, *) '         material property cards in referenc
     &e IN file.'
          write (istdo, *) 'The excess entries will not be offerred as d
     &efault values.'
          write (istdo, *) '********************************************
     &********************'
          numat = maxmat
        else
        end if
c
        read (irefu, *) (visc(k), k = 1, numat)
        read (irefu, *) (alam(k), k = 1, numat)
        read (irefu, *) (diff(k), k = 1, numat)
        if ((icode .eq. iddc) .or. (icode .eq. idds)) then
          read (irefu, *) (diffb (k), k = 1, numat)
        else
        end if
        read (irefu, *) (ra(k), k = 1, numat)
        if ((icode .eq. iddc) .or. (icode .eq. idds) .or. 
     &      (icode .eq. ichn) .or. (icode .eq. idef)) then
          read (irefu, *) (rab(k), k = 1, numat)
        else
        end if
        read (irefu, *) (dmhu(k), k = 1, numat)
        read (irefu, *) (estar(k), k = 1, numat)
        read (irefu, *) (toff(k), k = 1, numat)
        if (icode .ne. ichn) then
          read (irefu, *) (vstar(k), k = 1, numat)
          read (irefu, *) (x2ref(k), k = 1, numat)
          read (irefu, *) (sigref(k), k = 1, numat)
          read (irefu, *) (viscut(k), k = 1, numat)
        else
        end if
      else
c
c... offer defaults given in 'deck.defaults.h'
c
      end if
c
200   write (istdo, *)
      numat = getint  ('Enter number of material property groups:', 
     &                  numat)
      if (numat .gt. maxmat)  then
        write (isdto, *)
        write (istdo, *) '**********************************************
     &******************'
        write (istdo, *) 'Error: Only ', maxmat, ' material properties s
     &upported!'
        write (istdo, *) '**********************************************
     &******************'
        go to 200
      else
      end if
c
c... surface flux / force
c
      if (iref) then
c
        if (numsuf .gt. maxsuf) then
          write (istdo, *)
          write (istdo, *) '********************************************
     &********************'
          write (isdto, *) 'Warning: More than the (arbitrary) limit of
     &maxsuf = ', maxsuf
          write (isdto, *) '         surface force / flux cards in refer
     &ence IN file.'
          write (istdo, *) 'The excess entries will not be offerred as d
     &efault values.'
          write (istdo, *) '********************************************
     &********************'
          numsuf = maxsuf
        else
        end if
c
        do 220 k = 1, numsuf
          read (irefu, *) nel(k), iside(k), fnorm(k), ftan(k), flux(k)
220     continue
      else
c
c... offer defaults given in 'deck.defaults.h'
c
      end if
c
240   write (istdo, *)
      numsuf = getint ('Enter number of surface force/flux cards:', 
     &                  numsuf)
      if (numsuf .gt. maxsuf)  then
        write (isdto, *)
        write (istdo, *) '**********************************************
     &******************'
        write (istdo, *) 'Error: Only ', maxsuf, ' surface flux / force 
     &cards supported!'
        write (istdo, *) '**********************************************
     &******************'
        go to 240
      else
      end if
c
      write (iunit, "(i4,i7,10i4)" ) ntype, numel, nen, nenl, numat, 
     &                               nedof, numsuf, nitp, implv, implt
c
      do 300 k = 1, numat
        write (istdo, *)
        write (istdo, *) '**********************************************
     &******************'
        write (istdo, *) 'Order materials such that higher number'
        write (isdto, *) 'materials may overwrite lower number materials
     &.'
        write (istdo, *) '**********************************************
     &******************'
        write (istdo, *)
        if (numat .gt. 1) then
          write (istdo, "('***   Properties of material group', i3, 
     &                '   ***')") k
        else
        end if
        visc(k)    = gtreal ('Enter Viscosity:                    ', 
     &                        visc(k))
        alam(k)    = gtreal ('Enter Penalty Parameter:            ', 
     &                        alam(k))
        diff(k)    = gtreal ('Enter Thermal Diffusivity:          ', 
     &                        diff(k))
        if ((icode .eq. iddc) .or. (icode .eq. idds)) then
          diffb(k) = gtreal ('Enter Compositional Diffusivity:    ', 
     &                        diffb(k))
        else
        end if
        ra(k)      = gtreal ('Enter Thermal Rayleigh Number:      ', 
     &                        ra(k))
        if ((icode .eq. iddc) .or. (icode .eq. idds) .or. 
     &      (icode .eq. ichn) .or. (icode .eq. idef)) then
          rab(k)   = gtreal ('Enter Compositional Rayleigh Number:', 
     &                        rab(k))
        else
        end if
        dmhu(k)    = gtreal ('Enter Internal Heating Number:      ', 
     &                        dmhu(k))
        if (ntimvs .eq. 1) then
          estar(k)   = gtreal ('Enter Activation Energy:            ', 
     &                          estar(k))
          toff(k)    = gtreal ('Enter Temp Law Offset:              ', 
     &                          toff(k))
          if ((icode .eq. istd) .or. (icode .eq. iddc)) then
            vstar(k)   = gtreal ('Enter Activation Volume:            ',
     &                          vstar(k))
            x2ref(k)   = gtreal ('Enter Reference x2:                 ',
     &                          x2ref(k))
            if (nnnit .ne. 1) then
              sigref(k) = gtreal ('Enter Reference Stress:             '
     &,                            sigref(k))
            else
            end if
            viscut(k) = gtreal ('Enter Viscosity Cut-Off:            ',
     &                           viscut(k))
          else
          end if
        else
        end if
300   continue
c
      write (iunit, 9300) (visc(k), k = 1, numat)
      write (iunit, 9300) (alam(k), k = 1, numat)
      write (iunit, 9300) (diff(k), k = 1, numat)
      if ((icode .eq. iddc) .or. (icode .eq. idds)) then
        write (iunit, 9300) (diffb(k), k = 1, numat)
      else
      end if
      if(addcomment)write (iunit,*) 'rayleigh number'
      write (iunit, 9300) (ra(k), k = 1, numat)
      if ((icode .eq. iddc) .or. (icode .eq. idds) .or.
     &    (icode .eq. ichn) .or. (icode .eq. idef)) then
        write (iunit, 9300) (rab(k), k = 1, numat)
      else
      end if
      write (iunit, 9300) (dmhu(k), k = 1, numat)
      write (iunit, 9300) (estar(k), k = 1, numat)
      write (iunit, 9300) (toff(k), k = 1, numat)
      if (icode .ne. ichn) then
        write (iunit, 9300) (vstar(k), k = 1, numat)
        write (iunit, 9300) (x2ref(k), k = 1, numat)
        write (iunit, 9300) (sigref(k), k = 1, numat)
        write (iunit, 9300) (viscut(k), k = 1, numat)
      else
      end if
c
c... Surface Force / Flux Cards
c
      do 400 k = 1, numsuf
        write (istdo, *)
        write (istdo, *) '***   Surface Force / Flux Card #', k, '   ***
     &'
        nel(k) = getint ('Enter element number to apply surface force or
     & flux:', nel(k))
        write (istdo, *) '1 => bottom'
        write (istdo, *) '2 => right'
        write (istdo, *) '3 => top'
        write (istdo, *) '4 => left'
        iside(k) = getint ('Enter side to apply surface force or flux:
     &       ', iside(k))
        fnorm(k) = gtreal ('Enter normal surface force:        
     &       ', fnorm(k))
        ftan(k)  = gtreal ('Enter tangential surface force:        
     &       ', ftan(k))
        flux(k)  = gtreal ('Enter heat flux:        
     &       ', flux(k))
c
        write (iunit, "(2(i6,1x),3(1pe12.5))") nel(k), iside(k), 
     &                                         fnorm(k), ftan(k), 
     &                                         flux(k)
400   continue
c
      close (iunit)
      close (irefu)
c
c... build geometry file
c
c... Coordinate Group Cards
c
      write (istdo, *)
      write (istdo, *) '************************************************
     &****************'
      write (istdo, *) 'Coordinate Group Cards'
      write (istdo, *) 'Not (yet) based on existing input deck or includ
     &e file settings!'
      write (istdo, *) '************************************************
     &****************'
c
      if (gref) then
c
c... skip coordinate group cards
c
        call skpsec (grefu, 4)
      else
      end if
c
      write (istdo, *)
      if (yes ('Is this grid irregular', 'n')) then
        isubg = getint ('Enter number of grid sub-groups:        ', 2)
        node1 = ncorbl
      else
        node1 = 1
        nxg   = nelx
        nzg   = nelz
        isubg = 1
      end if

      do 500 i = 1, isubg
        if (isubg .gt. 1) then
          write (istdo, *)
          write (istdo, *) '***   Sub Group #', i, '   ***'
          nxg  = getint('Enter number of x1 elements in subgroup:',nelx)
          nzg  = getint('Enter number of x2 elements in subgroup:',nelz)
          node1= getint('Enter node number of lower left corner: ',node1
     &)
          xmin = gtreal('Enter x1 minimum value:   ', xmin)
          xmax = gtreal('Enter x1 maximum value:   ', xmax)
          zmin = gtreal('Enter x2 minimum value:   ', zmin)
          zmax = gtreal('Enter x2 maximum value:   ', zmax)
        else
        end if
c
        node2   = (nxg * incx) + node1
        node3   = (nxg * incx) + (nzg * incz) + node1
        node4   = (nzg * incz) + node1
c
        write (gunit, 9500)  node1, 4     , xmin , zmin 
        write (gunit, 9500)  node2, 1     , xmax , zmin 
        write (gunit, 9500)  node3, 1     , xmax , zmax
        write (gunit, 9500)  node4, 1     , xmin , zmax
        write (gunit, 9550) nxg, incx, nzg, incz
c
        if (lmovie) then
          write (munit, 9500)  node1, 4     , xmin , zmin 
          write (munit, 9500)  node2, 1     , xmax , zmin 
          write (munit, 9500)  node3, 1     , xmax , zmax
          write (munit, 9500)  node4, 1     , xmin , zmax
          write (munit, 9550) nxg, incx, nzg, incz
        else
        end if
500   continue
c
      write (gunit, 9500) 0, 0, 0.0, 0.0
      if (lmovie) then
        write (munit, 9500) 0, 0, 0.0, 0.0
      else
      end if
c
c... Velocity Boundary Condition Group - only simplest option 
c      implemented. Not (yet) based on existing input deck or include 
c      file settings!
c
      if (gref) then
c
c... skip velocity boundary condition group cards
c
        call skpsec (grefu, 4)
      else
      end if
c
c... explicitly pin velocities to facilitate any subsequent user changes
c
      if ((ivxb .eq. ifix) .or. (ivzb .eq. ifix)) then
        write (gunit, 9500) ncorbl, 2, vxbotl, vzbotl
        write (gunit, 9500) ncorbr, 0, vxbotr, vzbotr
        write (gunit, 9550) nelx, incx, 0, 0
      else
      end if       
c
      if ((ivxt .eq. ifix) .or. (ivzt .eq. ifix)) then
        write (gunit, 9500) ncortl, 2, vxtopl, vztopl
        write (gunit, 9500) ncortr, 0, vxtopr, vztopr
        write (gunit, 9550) nelx, incx, 0, 0
      else
      end if       
c
      write (gunit, 9500) 0, 0, 0.0, 0.0
c
c... Temperature/Composition Boundary Condition Group
c
      write (istdo, *)
      write (istdo, *) '************************************************
     &****************'
      write (istdo, *) 'Temperature/Composition Boundary Condition Cards
     &'
      write (istdo, *) 'Not (yet) based on existing input deck or includ
     &e file settings!'
      write (istdo, *) '************************************************
     &****************'
c
      do 600 loop = 1, nfield
c
c... make bottom and top fixed temperature/composition boundary 
c      conditions
c
c... bottom edge
c
        write (istdo, *)
c
        if (loop .eq. 1) then
          bot = gtreal ('Enter bottom temperature value:', tbot)
        else
          bot = gtreal ('Enter bottom composition value:', bbot)
        end if
c
        write (gunit, 9500) ncorbl, 2, bot
        write (gunit, 9500) ncorbr, 0, bot
        write (gunit, 9550) nelx, incx, 0, 0
c
c... top edge
c
        if (loop .eq. 1) then
          top = gtreal ('Enter top temperature value:   ', ttop)
        else
          top = gtreal ('Enter top composition value:   ', btop)
        end if
c
        write (gunit, 9500) ncortl, 2, top
        write (gunit, 9500) ncortr, 0, top
        write (gunit, 9550) nelx, incx, 0, 0
c
        write (gunit, 9500) 0, 0, 0.0
600   continue
c
c... Element Connectivity Generation Group Cards
c
      if (numat .gt. 1) then
        write (istdo, *)
        write (istdo, *) '**********************************************
     &******************'
        write (istdo, *) 'Element Connectivity Generation Group Cards'
        write (istdo, *) 'Not (yet) based on existing input deck or incl
     &ude file settings!'
        write (istdo, *) '**********************************************
     &******************'
      else
      end if
c
      do 750 matno = 1, numat
        if (numat .eq. 1) then
c
          ielnu = 1
          node1 = 1
          node2 = node1 + incx
          node3 = node2 + incz
          node4 = node1 + incz
          nxg = nelx
          nzg = nelz
c
          write (gunit, 9550) ielnu, 1, matno, node1, node2, node3,
     &                        node4
          write (gunit, 9550) nxg, incx-1, incx, nzg, incz, incz
c
          if (lmovie) then
            write (munit, 9550) ielnu, 1, matno, node1, node2, node3,
     &                          node4
            write (munit, 9550) nxg, incx-1, incx, nzg, incz, incz
          else
          end if
c
        else
c
          dfault = 'y'
          done = .false.
c
          do 700 while (.not. done)
c
            write (istdo, *)
            write (istdo, *) '***   Material Group #', matno, '   ***'
            node1 = getint('Enter lower left node of material group:    
     & '                  , node1)
            nxg   = getint('Enter number of elements in the x1 direction
     &:'                  , nxg)
            nzg   = getint('Enter number of elements in the x2 direction
     &:'                  , nzg)
            dum   = float (node1 / incx)
            irow  = int (dum)
            icol  = mod (node1, incx)
            ielnu = irow * (incx-1) + (icol * incz)
            node2 = node1 + incx
            node3 = node2 + incz
            node4 = node1 + incz
c
            write (gunit, 9550) ielnu, 1, matno, node1, node2, node3,
     &                          node4
            write (gunit, 9550) nxg, incx-1, incx, nzg, incz, incz
c
            if (lmovie) then
              write (munit, 9550) ielnu, 1, matno, node1, node2, node3,
     &                            node4
              write (munit, 9550) nxg, incx-1, incx, nzg, incz, incz
            else
            end if
c
            write (istdo, *)
            if (yes ('Finished with this material', dfault)) then
              done = .true.
            else
              done = .false.
            end if
c
700       continue
c
        end if
c
750   continue
c
      write (gunit, 9550) 0, 0, 0, 0, 0, 0, 0
      if (lmovie) then
        write (munit, 9550) 0, 0, 0, 0, 0, 0, 0
      else
      end if
c
      write(gunit,*)
      close (gunit)
      close (grefu)
      if (lmovie) then
        close (munit)
      else
      end if
c
      write (runit, 9800) ifile (1:kblnk(ifile))
      write (runit, 9800) gfile (1:kblnk(gfile))
      write (runit, 9800) outfn (1:kblnk(outfn)) // sep // 
     &                   suffix (1:kblnk(suffix))
      write (runit, 9800) rinfn (1:kblnk(rinfn)) // sep //
     &                   suffix (1:kblnk(suffix))
      write (runit, 9800) routfn (1:kblnk(routfn)) // sep //
     &                   suffix (1:kblnk(suffix))
      write (runit, 9800) tsfn (1:kblnk(tsfn)) // sep //
     &                   suffix (1:kblnk(suffix))
      write (runit, 9800) fldfn (1:kblnk(fldfn)) // sep //
     &                   suffix (1:kblnk(suffix))
      write (runit, 9800) meanfn (1:kblnk(meanfn)) // sep //
     &                   suffix (1:kblnk(suffix))
      write (runit, 9800) strfn (1:kblnk(strfn)) // sep //
     &                   suffix (1:kblnk(suffix))
      write (runit, 9800) cordfn (1:kblnk(cordfn)) // sep //
     &                   suffix (1:kblnk(suffix))
c
      if (icode .eq. ichn) then
        write (runit, 9800) richfn (1:kblnk(fldfn)) // sep //
     &                     suffix (1:kblnk(suffix))
        write (runit, 9800) rochfn (1:kblnk(fldfn)) // sep //
     &                     suffix (1:kblnk(suffix))
        write (runit, 9800) chnfn (1:kblnk(fldfn)) // sep //
     &                     suffix (1:kblnk(suffix))
      else
      end if
c
      close (runit)
c
9000  format (i7, 2i4, 2i5, 6i4, i6, 4i4, i9, 2i4)
9025  format (3i6, 1pe12.5)
9050  format (3 (i6, 1x), 5x, 3 (i3, 1x))
9300  format (10 (1pe12.5, 1x))
9500  format (2 (i6, 1x), 5x, 3 (f12.6, 1x))
9550  format (11 (i6, 1x))
9800  format (a)
c
      stop
      end
