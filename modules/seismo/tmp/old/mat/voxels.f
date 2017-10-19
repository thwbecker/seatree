	subroutine path_3d(n1layer,nsqrs,nsqtot,nlatzones,
     &ddep,ddel,EPLAT,EPLON,STLAT,STLON,timeisc,
     &rboundary,rescale,nlay,NCOLTOTAL,nto2,nto3,
     &nto4,nto5,nto8,nto9,cutoff,isotro,eq_incr,nlatzomax)
	implicit real*8 (a-h,o-z)
	real*4 eq_incr
c	parameter(ndimmax1=300,ndimmax2=3000)
	parameter(ndimmax1=500,ndimmax2=5000)
	common/lapo/raddown(ndimmax2),radup(ndimmax2),kdown,kup,time(ndimmax2),
     &dedown(ndimmax2),deup(ndimmax2),culdel(ndimmax2),culray(ndimmax2),iswtc
	common/tempo/premtime
c--cmb structure kernel (Morelli & Dziewonski, 1987)
	COMMON/CMBKER/cmbkernl
c--integrated anisotropic kernels--vph,vpv,vsh,vsv,eta, respectively 
	common/lapoanis/dqdvph(ndimmax2),dqdvpv(ndimmax2),dqdvsh(ndimmax2),
     &dqdvsv(ndimmax2),dqdeta(ndimmax2)
	real*4 deltami,deltama
	common/deltastats/deltami,deltama

	DIMENSION PHT(ndimmax1),IDX(ndimmax1)
	DIMENSION jsqrm(ndimmax1),ddelm(ndimmax1),radii2(ndimmax1)
	dimension indgrd(ndimmax1),indgrd2(ndimmax1),idxxx(1000),idxx(1000)
	double precision drsple
	real*4 rtarray(ndimmax1)
	real*4 rtarray2(ndimmax1)
	dimension numcell(ndimmax1)
	real*4 wei
	real*4 eplat,eplon,stlat,stlon
	real*8 phtime(ndimmax1)
	real*8 phdqdvph(ndimmax1),phdqdvpv(ndimmax1)
	real*8 dedown2(ndimmax2),deup2(ndimmax2),raddown2(ndimmax2),radup2(ndimmax2)
	real*8 coeffup(3,ndimmax2),coeffdown(3,ndimmax2),F(3,ndimmax2)
	real*8 rboundary(0:20),radii(ndimmax2),radintr(ndimmax2),delintr(ndimmax2)
	real*8 rrr,ddd,coeff(3,ndimmax2),f2(3,ndimmax2)
	real*8 radup3(ndimmax2),raddown3(ndimmax2)
	real*8 coefft(3,ndimmax2)
	real*8 coeffani1(3,ndimmax2),coeffani2(3,ndimmax2)
	real*4 rescale(20),cutoff
	character*20 outexc
	real*4 outwri
c	integer nsqrs(nlatzones),nsqtot(nlatzones+1)
	integer nsqrs(nlatzomax),nsqtot(nlatzomax+1)

c      DATA NSQRS/4,10,16,22,28,34,38,44,48,54,58,60,64,66,68,70,72,72,
c     .           72,72,70,68,66,64,60,58,54,48,44,38,34,28,22,16,10,4/

	parameter(pi=3.14159265358979d0)
	parameter(radian=180.d0/pi)
	DATA GEOCO/0.993277/

c	PARAMETER(N1LAYER=1656)
	NALLAYERS=NLAY*N1LAYER

cTEST
c        write(99,*)">",ddel

	nthrown=0
	iscrup=0
c--number of entries of the kernel array
	nkern=kup+iswtc
c--all the arrays given by ttpath (common block "lapo")
c--need to be reorganized to be used by this subroutine:
	do iti=kdown+kup+1,2,-1
	   time(iti)=time(iti-1)
	   dqdvph(iti)=dqdvph(iti-1)
	   dqdvpv(iti)=dqdvpv(iti-1)
	   dqdvsh(iti)=dqdvsh(iti-1)
	   dqdvsv(iti)=dqdvsv(iti-1)
	   dqdeta(iti)=dqdeta(iti-1)
cTEST
c        write(99,*)time(iti),dqdvph(iti),dqdvpv(iti)

	enddo
	time(1)=0.
	dqdvph(1)=0.
	dqdvpv(1)=0.
	dqdvsh(1)=0.
	dqdvsv(1)=0.
	dqdeta(1)=0.

	time(kdown+kup+2)=premtime

	rsource=6371.d0-ddep
	do irdw=kdown+1,2,-1
	   raddown(irdw)=raddown(irdw-1)
	   dedown(irdw)=dedown(irdw-1)
	enddo
	raddown(1)=rsource
	dedown(1)=0.
	kdown=kdown+1

	radup(kup+1+iswtc)=6371.d0
	deup(kup+1+iswtc)=ddel*pi/180.d0
	kup=kup+1

	do irdw=1,kdown
	   culdel(irdw)=dedown(irdw)
	   culray(irdw)=raddown(irdw)
	enddo

	do irup=1,kup
	   culdel(irup+kdown)=deup(irup+iswtc)
	   culray(irup+kdown)=radup(irup+iswtc)
	enddo

	DO KW=1,KdowN+kup
	   if((culray(kw).le.0.).or.(culdel(kw).lt.0.))then
c	      write(94,*)'screwed up path'
	      iscrup=1
	      nto3=nto3+1
	      goto4366
	   endif
	enddo
	rbottom=radup(1)
	dbottom=deup(1)

	DTH=PI/float(nlatzones)
c--dth is 5 degrees: increment in latitude.
	PI2=PI*2.
	FORPI=4.*PI
	nth1=nlatzones+1
c--convert to radians & correct latitude for geocentric coordinates
       EPLgc=radian*ATAN(GEOCO*TAN(dble(EPLAT)/radian))
      TH1=(90.-dble(EPLgc))/RADIAN
      PH1=(dble(EPLON))/RADIAN
       stLgc=radian*ATAN(GEOCO*TAN(dble(stlat)/radian))
      TH2=(90.-dble(STLgc))/RADIAN
      PH2=(dble(STLON))/RADIAN
      STH1=dSIN(TH1)
      STH2=dSIN(TH2)
      CTH1=dCOS(TH1)
      CTH2=dCOS(TH2)
      SPH1=dSIN(PH1)
      SPH2=dSIN(PH2)
      CPH1=dCOS(PH1)
      CPH2=dCOS(PH2)
c--find out the coordinates of the "north pole of the ray" knowing
c--the coordinates of the source and of the receiver:
      CPH21= CPH1*CPH2+SPH1*SPH2
      SPH21=SPH2*CPH1-SPH1*CPH2
      CDEL=STH1*STH2*CPH21+CTH1*CTH2
	if(cdel.eq.-1.)then
	   print*,'delta too large; equations break down'
	   nto9=nto9+1
	   goto4366
	endif
      CCAPTH=STH1*STH2*SPH21/dSQRT(1.-CDEL*CDEL)
      SCAPTH=dSQRT(1.-CCAPTH*CCAPTH)
      capth=datan2(scapth,ccapth)
      SCAPPH=CTH1*STH2*CPH2-CTH2*STH1*CPH1
      CCAPPH=STH1*CTH2*SPH1-STH2*CTH1*SPH2
      CAPPH=dATAN2(SCAPPH,CCAPPH)
      SCAPPH=dSIN(CAPPH)
      CCAPPH=dCOS(CAPPH)
c--capth (capital theta) and capph (capital phi) are now colatitude and
c--longitude of the "north pole" of the ray's great circle w/ respect to
c--the true north pole. cdel is cosine of epi. distance.

c--determine del (epicentral distance) from cosine(del)
	DEL=dATAN2(dSQRT(1.-CDEL*CDEL),CDEL)
c--this delta is in general slightly different from the one determined by
c--ttmain; in fact here an exact formula is used while ttmain works by
c--successive approximations.
	CPHSP=CCAPTH*STH1*(CPH1*CCAPPH+SPH1*SCAPPH)-SCAPTH*CTH1
	SPHSP=STH1*(SPH1*CCAPPH-CPH1*SCAPPH)
	PHSP=dATAN2(SPHSP,CPHSP)
c--phsp is (?) probably the longitude of the source ON the great circle;
c--how is it defined = what is its zero?
c--phsp is not used to compute the intersections with the grid-
c--those are determined with respect to an arbitrary (?) point;
c--it's used later when jsqre is figured out

c--the following determines the max & min latitudes reached 
c--by the great circle
      thet=capth
      if(capth.gt.0.5*pi) thet=pi-capth
      thmin=0.5*pi-thet
      thmax=0.5*pi+thet

c--INITIALIZE
      lat_zone=0
      IENT=0
      IF(SCAPTH.EQ.0) GOTO 10

	do kp=1,ndimmax1
	   pht(kp)=0.
           indgrd(kp)=0
	   idx(kp)=0
	   radii(kp)=0.
	   radii2(kp)=0.
	enddo

cTEST	
c	open(888,file="merid.int")
c	open(889,file="paral.int")
c	print*,"saving latitudinal and longitudinal intersections"

c--FIND deltas of INTERSECTIONS WITH THE horizontal GRID:
C--(they depend only on source and station locations)
	DO 20 I=2,nth1
         lat_zone=lat_zone+1
c--th is cumulative latitudinal increment
         TH=dble(FLOAT(I-1))*DTH
         CTH=dCOS(TH)
c--cumulative delta(pht): cos (only spherical trigonometry eq.)
         CPHT=-CTH/SCAPTH
         CPHT2=CPHT*CPHT
c--skip latitude zone if it's out of the lat. range of the ray:
         if(th.lt.thmin.or.th-dth.gt.thmax) go to 20
c--the following takes care of the last latitude zone, for which
c--there is no intersection with the "th" parallel, only with "th-dth":
         IF(CPHT2.GT.1.)then
	    GOTO 21
	 endif
c--if everything's ok, store in pht the cumulative delta
c--these are "latitudinal intersections"
c--in general a ray crosses a parallel in 2 points 
         IENT=1+IENT
         PHT(IENT)=dATAN2(dSQRT(1.-CPHT2),CPHT)
         indgrd(ient)=1
         spht=dsin(pht(ient))
         IENT=1+IENT
         PHT(IENT)=-PHT(IENT-1)
cTEST--------------------------
c	pht1=pht(ient)
c      CPHT_1=dCOS(PHT1)
c      SPHT_1=DSIN(PHT1)
c      CTH_1=-CPHT_1*SCAPTH
c      TH_1=dATAN2(dSQRT(1.-CTH_1*CTH_1),CTH_1)
c      CPH_1=CPHT_1*CCAPTH*CCAPPH-SPHT_1*SCAPPH
c      SPH_1=CPHT_1*CCAPTH*SCAPPH+SPHT_1*CCAPPH
c      IF(SPH_1.EQ.0..AND.CPH_1.EQ.0.)GOTO 2873
c        PH_1=dATAN2(SPH_1,CPH_1)
c      if(ph_1.lt.0.) ph_1=ph_1+pi2
c      if(ph_1.gt.pi2) ph_1=ph_1-pi2
c      GOTO 2874
c 2873 PH=0.
c 2874 continue
c	write(889,*)ph_1*180/pi,90.-(th_1*180/pi),ient,indgrd(ient)
c------------------------------
cTEST--------------------------
c	pht1=pht(ient-1)
c      CPHT_1=dCOS(PHT1)
c      SPHT_1=DSIN(PHT1)
c      CTH_1=-CPHT_1*SCAPTH
c      TH_1=dATAN2(dSQRT(1.-CTH_1*CTH_1),CTH_1)
c      CPH_1=CPHT_1*CCAPTH*CCAPPH-SPHT_1*SCAPPH
c      SPH_1=CPHT_1*CCAPTH*SCAPPH+SPHT_1*CCAPPH
c      IF(SPH_1.EQ.0..AND.CPH_1.EQ.0.)GOTO 3873
c        PH_1=dATAN2(SPH_1,CPH_1)
c      if(ph_1.lt.0.) ph_1=ph_1+pi2
c      if(ph_1.gt.pi2) ph_1=ph_1-pi2
c      GOTO 3874
c 3873 PH=0.
c 3874 continue
c	write(889,*)ph_1*180/pi,90.-(th_1*180/pi),ient-1,indgrd(ient-1)
c------------------------------

  21     numlong=nsqrs(lat_zone)
         dphi=pi2/dble(float(numlong))

c--dphi will be longitude increment
c--since the longitude increment changes depending on the lat.zone,
c--we have to nest the following loop, that computes the intersections
c--with meridians, inside the loop over latitude zones...

c--now in each latitudinal zone find the epicentral distances
c--of intersections of the great circle with the longitudinal
c--boundaries of the grid.
          DO 40 j=1,numlong
c--ph is cumulative longitudinal increment
            ph=dble(float(j-1))*dphi
            angr=ph-capph
            thlo=atan(-ccapth/(scapth*cos(angr)))
            if(thlo.lt.0.)thlo=pi+thlo
            if(thlo.gt.th-dth.and.thlo.lt.thmin)then
	       goto40
	    endif
            if(thlo.lt.th+dth.and.thlo.gt.thmax)then
	       goto40
	    endif
            if(thlo.gt.th.or.thlo.lt.th-dth)then
	       goto40
	    endif
            SPH=dSIN(PH)
            CPH=dCOS(PH)
            IENT=IENT+1
            PHT(IENT)=dATAN2(CCAPTH*(SPH*CCAPPH-CPH*SCAPPH),
     &                CPH*CCAPPH+SPH*SCAPPH)
        indgrd(ient)=2
            IF(PHT(IENT).GT.PI) then
	       PHT(IENT)=PHT(IENT)-PI2
	    endif
cTEST--------------------------
c	pht1=pht(ient)
c      CPHT_1=dCOS(PHT1)
c      SPHT_1=DSIN(PHT1)
c      CTH_1=-CPHT_1*SCAPTH
c      TH_1=dATAN2(dSQRT(1.-CTH_1*CTH_1),CTH_1)
c      CPH_1=CPHT_1*CCAPTH*CCAPPH-SPHT_1*SCAPPH
c      SPH_1=CPHT_1*CCAPTH*SCAPPH+SPHT_1*CCAPPH
c      IF(SPH_1.EQ.0..AND.CPH_1.EQ.0.)GOTO 1873
c        PH_1=dATAN2(SPH_1,CPH_1)
c      if(ph_1.lt.0.) ph_1=ph_1+pi2
c      if(ph_1.gt.pi2) ph_1=ph_1-pi2
c      GOTO 1874
c 1873 PH=0.
c 1874 continue
c	write(888,*)ph_1*180/pi,90.-(th_1*180/pi),ient,indgrd(ient)
c------------------------------

   40 CONTINUE
   20 CONTINUE
cTEST
c	close(888)
c	close(889)

   10 CONTINUE

c--now the CUMULATIVE deltas corresponding to each square crossed by
c--the ray have been computed, and stored in pht.

      DO 60 I=1,IENT
	PHT(I)=dMOD(PHT(I)-PHSP+FORPI,PI2)
60    continue
	del=ddel*pi/180.d0
      IENT=IENT+1
      PHT(IENT)=DEL
      indgrd(ient)=4
      IENT=IENT+1
      PHT(IENT)=0.
      indgrd(ient)=4
c--at this point pht = cumulative delta (in an arbitrary order)

c--DETERMINE THE RADII OF THE INTERSECTIONS WITH THE HORIZONTAL GRID
c-------------expand radius(delta) over splines.
	call drspln(1,kdown+kup,culdel,culray,coeff,f2)

c--interpolate the function radius(delta) at the deltas already stored in pht
c--notice we make use of the fact that the source is the last element stored
c--in pht.
	isou=0
	ingridold=0
	do ip=1,ndimmax1
	   if((pht(ip).le.del).and.(isou.eq.0))then
	      ddd=pht(ip)
	      radii(ip)=drsple(1,kup+kdown,culdel,culray,coeff,ddd)
	      if((radii(ip).gt.6371.).or.(radii(ip).lt.0.))then
	       iscrup=1
c	       write(94,*)'screwed up interpolation',radii(ip)
	       nto2=nto2+1
	       goto4366
	      endif
ccccc	      iiiii=ieee_flags("get","exception","invalid",outexc)
ccccc	      if(outexc.eq.'invalid')then
ccccc	       iscrup=1
c	       write(94,*)'screwed up interpolation',radii(ip),testnan
ccccc	       nto2=nto2+1
ccccc	       iiiii=ieee_flags("clear","exception","all",out)
ccccc	       goto4366
ccccc	      endif
	      if((indgrd(ip).eq.4).and.(ingridold.eq.4))isou=1
	      ingridold=indgrd(ip)
	   endif
	enddo

c----------------------what is this for?
	ktota=kdown+kup
	do intr=2,kdown+kup
	   ddinc=(culdel(intr)-culdel(intr-1))/2.
	   ddd=culdel(intr-1)+ddinc
	   radintr(intr-1)=drsple(1,kup+kdown,culdel,culray,coeff,ddd)
	   delintr(intr-1)=ddd
	enddo

c--to make "radial" and "grid" interpolation consistent with each other,
c--add to dedown, deup, raddown, radup the interpolated points corresponding
c--to "grid" intersections (pht).
	isou=0
	do inte=1,ndimmax1
	   if((pht(inte).lt.del).and.(isou.eq.0))then
	      if(pht(inte).lt.deup(1))then
	         kdown=kdown+1
	         dedown(kdown)=pht(inte)
	         raddown(kdown)=radii(inte)
	      elseif(pht(inte).gt.deup(1))then
	         kup=kup+1
	         deup(kup+iswtc)=pht(inte)
	         radup(kup+iswtc)=radii(inte)
	      endif
	   endif
	   if(pht(inte).eq.0.)isou=1
	enddo

c--repeat the same procedure for the arbitrarily added points (radintr, delintr)
	do inte=1,ktota-1
	      if(delintr(inte).lt.deup(1))then
	         kdown=kdown+1
	         dedown(kdown)=delintr(inte)
	         raddown(kdown)=radintr(inte)
	      elseif(delintr(inte).gt.deup(1))then
	         kup=kup+1
	         deup(kup+iswtc)=delintr(inte)
	         radup(kup+iswtc)=radintr(inte)
	      endif
	enddo

c---------------------------------------sort
	CALL RSOINC(dedown,kdown,idxx)

	do ki=1,kdown
	   idd=idxx(ki)
	   raddown3(ki)=raddown(idd)
	enddo

	CALL RSOINC(deup,kup+iswtc,idxxx)

	do ki=1,kup+iswtc
	   idd=idxxx(ki)
	   radup3(ki)=radup(idd)
	enddo

c--clean up the function delta(radius) in order to interpolate it.
	kount=1
	raddown2(1)=raddown3(1)
	dedown2(1)=dedown(1)
	do im=2,kdown
	   if(dabs(raddown3(im)-raddown3(im-1)).gt.0.1d0)then
	      kount=kount+1
	      raddown2(kount)=raddown3(im)
	      dedown2(kount)=dedown(im)
	   elseif(dabs(raddown3(im)-raddown3(im-1)).le.0.1d0)then
	      kount=kount+1
	      raddown2(kount)=raddown2(kount-1)
	      dedown2(kount)=dedown(im)
	   endif
	enddo
	kdown=kount

	kount=1
	radup2(1)=radup3(1)
	deup2(1)=deup(1)
	do im=2,kup+iswtc
	   if(dabs(radup3(im)-radup3(im-1)).gt.0.1d0)then
	      kount=kount+1
	      radup2(kount)=radup3(im)
	      deup2(kount)=deup(im)
	   elseif(dabs(radup3(im)-radup3(im-1)).le.0.1d0)then
	      kount=kount+1
	      radup2(kount)=radup2(kount-1)
	      deup2(kount)=deup(im)
	   endif
	enddo
	kup=kount-iswtc

c--add to pht the (cumulative) deltas corresponding to intersections of
c--the ray with boundaries between layers.
c--it's done only for minor arc.
c--In order to do that, rayup and raydown must be interpolated.
c--first, expand them over splines.
	call drspln(1,kdown,raddown2,dedown2,coeffdown,f)
	call drspln(1,kup+iswtc,radup2,deup2,coeffup,f)

c--now we should have the function delta(radius), expressed as a cubic
c--splines expansion (coeffdown for the downgoing part and coeffup
c--for the upgoing).
c--now interpolate:
cTEST
c	open(734,file="radii.txt")

	kr=0
	do il=1,nlay
	   rrr=rboundary(il)
	   if((rrr.lt.rsource).and.(rrr.ge.rbottom))then
	      kr=kr+1
	      ient=ient+1
	      pht(ient)=drsple(1,kdown,raddown2,dedown2,coeffdown,rrr)
              indgrd(ient)=3
	      radii(ient)=rrr
cTEST
c	write(734,*)pht(ient),rrr,rbottom

	   endif
	enddo

	do il=1,nlay
	   rrr=rboundary(il)
	   if(rrr.ge.rbottom)then
	      ient=ient+1
	      kr=kr+1
	      pht(ient)=drsple(1,kup+iswtc,radup2,deup2,coeffup,rrr)
              indgrd(ient)=3
	      radii(ient)=rrr
cTEST
c	write(734,*)pht(ient),rrr,rbottom

	   endif
	enddo

cTEST
c	close(734)

c--at this point I have pht and radii and they are complete: sort
      CALL RSOINC(PHT,IENT,IDX)

c--rsoinc has now sorted pht by increasing cumulative delta.
c--each element of idx gives the position that the corresponding
c--element of pht had before sorting.
c------------------------sort radii
	do ki=1,ient
	   idd=idx(ki)
	   radii2(ki)=radii(idd)
           indgrd2(ki)=indgrd(idd)
	enddo

4567	continue

      PHT(IENT+1)=PI2
      ientm=0
      ientg=0
      cdelta=0.

	do ic=1,ient
	   if(pht(ic).eq.del)radrec=radii2(ic)
	enddo
	iscrew=0
	if(dabs(radrec-radup2(kup+iswtc)).gt..5)then
	   iscrew=1
	   nto5=nto5+1
	   goto4366
	endif

c--determine incremental travel time corresponding to array pht
	call drspln(1,ktota,culdel,time,coefft,f2)
	do ip=1,ndimmax1
	   if(pht(ip).le.del)then
	      ddd=pht(ip)
	      phtime(ip)=drsple(1,ktota,culdel,time,coefft,ddd)
	      if(pht(ip).eq.del)goto7889
	   endif
	enddo
7889	continue
c--determine incremental anisotropic kernel vph
	call drspln(1,ktota-1,culdel,dqdvph,coeffani1,f2)
	do ip=1,ient
	   if(pht(ip).le.del)then
	      ddd=pht(ip)
	      phdqdvph(ip)=drsple(1,ktota-1,culdel,dqdvph,coeffani1,ddd)
	      if(pht(ip).eq.del)goto7890
	   endif
	enddo
7890	continue
c--determine incremental anisotropic kernel vpv
	call drspln(1,ktota-1,culdel,dqdvpv,coeffani2,f2)
	do ip=1,ient
	   if(pht(ip).le.del)then
	      ddd=pht(ip)
	      phdqdvpv(ip)=drsple(1,ktota-1,culdel,dqdvpv,coeffani2,ddd)
	      if(pht(ip).eq.del)goto7891
	   endif
	enddo
7891	continue

cTEST
c        write(98,*)">",ddel
c        do ip=1,ient
c        write(98,*)phtime(ip),phdqdvph(ip),phdqdvpv(ip)
c        enddo
c        write(97,*)">",ddel
        coutker1=0.
        coutker2=0.

c----------------------determine indices of blocs crossed by ray-path
	NCOL=0
	idiscard=0

	indcmb1=0
	indcmb2=0
	itrans=0

c--loop over the entire great circle
      DO 50 I=1,IENT
      I1=1+I
c--calculate coordinates of mid-point of each ray segment
      PHTT=.5*(PHT(I)+PHT(I1))+PHSP
      CPHT=dCOS(PHTT)
      SPHT=DSIN(PHTT)
      CTH=-CPHT*SCAPTH
      TH=dATAN2(dSQRT(1.-CTH*CTH),CTH)
      CPH=CPHT*CCAPTH*CCAPPH-SPHT*SCAPPH
      SPH=CPHT*CCAPTH*SCAPPH+SPHT*CCAPPH
      IF(SPH.EQ.0..AND.CPH.EQ.0.)GOTO 9873
      PH=dATAN2(SPH,CPH)
      if(ph.lt.0.) ph=ph+pi2
      if(ph.gt.pi2) ph=ph-pi2
      GOTO 9874
 9873 PH=0.
 9874 CONTINUE
      XLAT=90.-TH*RADIAN
      XLON=PH*RADIAN
c--then the function isqre gives the "horizontal" index
	jsqre=isqre(xlat,xlon,nsqrs,nsqtot,nlatzones,n1layer,eq_incr)
c--incremental delta
	RD=(PHT(I1)-PHT(I))
c--incremental travel time
	rt=(phtime(i1)-phtime(i))
c--integrals of the anisotropic kernels in a block
	outker1=phdqdvph(i1)-phdqdvph(i)
	outker2=phdqdvpv(i1)-phdqdvpv(i)

c--if the segment belongs to the minor arc also find out the layer
      IF(PHT(I1).LE.DEL) then
	ientm=ientm+1
	jsqrm(ientm)=jsqre
	ddelm(ientm)=rd
	layenold=layen
	rady=(radii2(i)+radii2(i1))/2.
	do ila=1,nlay
	   if((rady.le.rboundary(ila-1)).and.(rady.ge.rboundary(ila)))then
	   layen=ila
	   endif
	enddo

cTEST
c	write(887,*)xlon,xlat,jsqre,layen

	if(rady.le.rboundary(nlay))layen=nlay+1
c--find the index of the CMB cell (if there is any) hit by the ray.
	IF((RADY.lt.RBOUNDARY(NLAY)).and.(itrans.eq.0))then
	   itrans=1
	   indcmb1=jsqre
	endif
	IF((RADY.gt.RBOUNDARY(NLAY)).and.(itrans.eq.1))then
	   itrans=0
	   indcmb2=jsqre
	endif
	IF(RADY.GE.RBOUNDARY(NLAY))then
c--the next "if" is a fix because there are repeated indexes.
         if( (ientm.gt.1.and.jsqrm(ientm).eq.jsqrm(ientm-1)) .and. 
     &   (layenold.eq.layen)) then
	    if((rd.ge.(0.001)).and.(rdold.ge.0.001))then
	       write(*,*)'REPEATED INDEX!', rdold
cTEST
c	stop

c--in this case the error is too big to be neglected - probably due
c--to inconsistencies between the radial and horizontal interpolation
c--so I'll have to throw away the datum
		idiscard=1
	        nto4=nto4+1
	        goto4366
	    endif
            ientm=ientm-1
c--in this case the error is small enough that we can save the datum anyway.
            ddelm(ientm)=ddelm(ientm)+rd
         endif
c--here in the anisotropic case things are different
	if(isotro.eq.2)then
c--count the number of non-zero elements in that ROW (row <--> path)
	NCOL=NCOL+2
C--STORE THE ANISOTROPIC KERNELS
c--remember: outker1 is the vph kernel, outker2 the vpv one.
	RTARRAY(NCOL-1)=OUTKER1
	RTARRAY(NCOL)=OUTKER2
c--store the cell (column) number:
	numcell(ncol-1)=((LAYEN-1)*n1layer)+JSQRE
	numcell(ncol)=  ((LAYEN-1)*n1layer)+JSQRE+nallayers
	elseif(isotro.eq.1)then
c--count the number of non-zero elements in that ROW (row <--> path)
	NCOL=NCOL+1
C--STORE THE ISOTROPIC KERNEL
c--remember: outker1 is the vph kernel, outker2 the vpv one.
	RTARRAY(NCOL)=OUTKER2+OUTKER1
	RTARRAY2(NCOL)=rt
cTEST
c	print*,outker1,outker2,outker1+outker2
c	if(outker1.gt.1..or.outker2.gt.1.)stop "positive value"
        coutker1=coutker1+outker1
        coutker2=coutker2+outker2
c        write(97,"(4(f11.5,2x))")phtime(i1),cdelta+rd,coutker1,coutker2


c--store the cell (column) number:
	numcell(ncol)=((LAYEN-1)*n1layer)+JSQRE
	endif
	endif
c1243	continue
	rdold=rd
	cdelta=cdelta+rd
	call extent(jsqre,XLAMIN,XLAMAX,XLOMIN,XLOMAX
     & ,nsqrs,nsqtot,nlatzones,eq_incr)
	endif
50	continue
c-----------------end of the loop over the entire great circle

cTEST
c	close(887)
c	print*,"loop regularly ended"
c	print*,idiscard,iscrew,iscrup
c	stop

c--store efficiently the matrix on disc
	if(((idiscard.eq.0).and.(iscrew.eq.0)).and.(iscrup.eq.0))then
c--assign a weight to this row
	wei=1.0
	DTME=abs(timeisc)
	if((dtme.gt.cutoff))then
	   wei=exp(-(dtme-cutoff))
	endif

c--statistics of epicentral distance
	if(ddel.lt.deltami)deltami=ddel
	if(ddel.gt.deltama)deltama=ddel

c--store nonzero elements+indices of a row of the rescaled matrix
	   do indice=1,ncol
	      ninte0=numcell(indice)
4141	      if(ninte0.gt.nallayers)then
	         ninte0=ninte0-nallayers
	      endif
	      if(ninte0.gt.nallayers)goto4141
	      ninte=((ninte0-1)/n1layer)
	      ninte=ninte+1
	      write(11,rec=ncoltotal+indice)
     &	      rtarray(indice)*wei
	      write(12,rec=ncoltotal+indice)numcell(indice)
cTEST
c	print*,indice,rtarray(indice),numcell(indice),wei
	
	   enddo
cTEST
c	pause

	if(cmbkernl.ne.0.)then
	   if((indcmb1.eq.0).and.(indcmb2.ne.0))then
	      print*,'error with the cmb stuff'
	      stop
	   elseif((indcmb1.ne.0).and.(indcmb2.eq.0))then
	      print*,'error with the cmb stuff'
	      stop
	   endif
c--want to write real*4 numbers
	   outwri=cmbkernl
	   write(11,rec=ncoltotal+ncol+1)outwri*wei
	   write(12,rec=ncoltotal+ncol+1)nallayers*isotro+indcmb1
	   write(11,rec=ncoltotal+ncol+2)outwri*wei
	   write(12,rec=ncoltotal+ncol+2)nallayers*isotro+indcmb2
	   ncol=ncol+2
	endif
	   ncoltotal=ncoltotal+ncol

c--store number of nonzero elements in this row
	   write(13,*)ncoltotal
c--rhs = measured total travel time - predicted total travel time
	   write(9,*)real(timeisc*wei)
	else
	   stop "problems writing matrix to file"
	endif

4366	continue
	return
	end

	function isqre(lat,lon,nsqrs,nsqtot,nlatzones,n,eq_incr)
c----finds the index of the square where (lat,lon) is
	real*8 lat,lon,loc_incr
	dimension nsqrs(nlatzones),nsqtot(nlatzones+1)
	lazone=(90.-lat)/eq_incr+1
	if((90.-lat).gt.180.)lazone=nlatzones
	if((90.-lat).gt.181.)stop "problems in function isqre"
	if(lazone.gt.nlatzones)then
	   print*,"problems in function isqre, latitude",lazone,lat
	   stop
	endif
	if(lon.lt.0.)lon=360.+lon
	loc_incr=360./float(nsqrs(lazone))
	isqre=(lon/loc_incr)+1
	isqre=isqre+nsqtot(lazone)
	if(isqre.gt.n)then
	   print*,"problems in function isqre, longitude"
	   stop
	endif
	RETURN
	END

C       FUNCTION DRSPLE
C$PROG DRSPLE
      DOUBLE PRECISION FUNCTION DRSPLE(I1,I2,X,Y,Q,S)
C
C C$C$C$C$C$ CALLS ONLY LIBRARY ROUTINES C$C$C$C$C$
C
C   RSPLE RETURNS THE VALUE OF THE FUNCTION Y(X) EVALUATED AT POINT S
C   USING THE CUBIC SPLINE COEFFICIENTS COMPUTED BY RSPLN AND SAVED IN
C   Q.  IF S IS OUTSIDE THE INTERVAL (X(I1),X(I2)) RSPLE EXTRAPOLATES
C   USING THE FIRST OR LAST INTERPOLATION POLYNOMIAL.  THE ARRAYS MUST
C   BE DIMENSIONED AT LEAST - X(I2), Y(I2), AND Q(3,I2).
C
C                                                     -RPB
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(1),Y(1),Q(3,1)
      DOUBLE PRECISION X,Y,Q,S
      DATA I/1/
      II=I2-1
C   GUARANTEE I WITHIN BOUNDS.
      I=MAX0(I,I1)
      I=MIN0(I,II)
C   SEE IF X IS INCREASING OR DECREASING.
      IF(X(I2)-X(I1))1,2,2
C   X IS DECREASING.  CHANGE I AS NECESSARY.
 1    IF(S-X(I))3,3,4
 4    I=I-1
      IF(I-I1)11,6,1
 3    IF(S-X(I+1))5,6,6
 5    I=I+1
      IF(I-II)3,6,7
C   X IS INCREASING.  CHANGE I AS NECESSARY.
 2    IF(S-X(I+1))8,8,9
 9    I=I+1
      IF(I-II)2,6,7
 8    IF(S-X(I))10,6,6
 10   I=I-1
      IF(I-I1)11,6,8
 7    I=II
      GO TO 6
 11   I=I1
C   CALCULATE RSPLE USING SPLINE COEFFICIENTS IN Y AND Q.
 6    H=S-X(I)
      DRSPLE=Y(I)+H*(Q(1,I)+H*(Q(2,I)+H*Q(3,I)))
      RETURN
      END


c================================================================
      SUBROUTINE DRSPLN(I1,I2,X,Y,Q,F)
C
C C$C$C$C$C$ CALLS ONLY LIBRARY ROUTINES C$C$C$C$C$
C
C   SUBROUTINE RSPLN COMPUTES CUBIC SPLINE INTERPOLATION COEFFICIENTS
C   FOR Y(X) BETWEEN GRID POINTS I1 AND I2 SAVING THEM IN Q.  THE
C   INTERPOLATION IS CONTINUOUS WITH CONTINUOUS FIRST AND SECOND
C   DERIVITIVES.  IT AGREES EXACTLY WITH Y AT GRID POINTS AND WITH THE
C   THREE POINT FIRST DERIVITIVES AT BOTH END POINTS (I1 AND I2).
C   X MUST BE MONOTONIC BUT IF TWO SUCCESSIVE VALUES OF X ARE EQUAL
C   A DISCONTINUITY IS ASSUMED AND SEPERATE INTERPOLATION IS DONE ON
C   EACH STRICTLY MONOTONIC SEGMENT.  THE ARRAYS MUST BE DIMENSIONED AT
C   LEAST - X(I2), Y(I2), Q(3,I2), AND F(3,I2).  F IS WORKING STORAGE
C   FOR RSPLN.
C                                                     -RPB
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION X,Y,Q,F
      DIMENSION X(1),Y(1),Q(3,1),F(3,1),YY(3)
      EQUIVALENCE (YY(1),Y0)
      DATA SMALL/1.D-10/,YY/0.D0,0.D0,0.D0/
      J1=I1+1
      Y0=0.D0
C   BAIL OUT IF THERE ARE LESS THAN TWO POINTS TOTAL.
      IF(I2-I1)13,17,8
 8    A0=X(J1-1)
C   SEARCH FOR DISCONTINUITIES.
      DO 3 I=J1,I2
      B0=A0
      A0=X(I)
      IF(DDABS((A0-B0)/DMAX1(A0,B0)).LT.SMALL) GO TO 4
 3    CONTINUE
 17   J1=J1-1
      J2=I2-2
      GO TO 5
 4    J1=J1-1
      J2=I-3
C   SEE IF THERE ARE ENOUGH POINTS TO INTERPOLATE (AT LEAST THREE).
 5    IF(J2+1-J1)9,10,11
C   ONLY TWO POINTS.  USE LINEAR INTERPOLATION.
 10   J2=J2+2
      Y0=(Y(J2)-Y(J1))/(X(J2)-X(J1))
      DO 15 J=1,3
      Q(J,J1)=YY(J)
 15   Q(J,J2)=YY(J)
      GO TO 12
C   MORE THAN TWO POINTS.  DO SPLINE INTERPOLATION.
 11   A0=0.
      H=X(J1+1)-X(J1)
      H2=X(J1+2)-X(J1)
      Y0=H*H2*(H2-H)
      H=H*H
      H2=H2*H2
C   CALCULATE DERIVITIVE AT NEAR END.
      B0=(Y(J1)*(H-H2)+Y(J1+1)*H2-Y(J1+2)*H)/Y0
      B1=B0
C   EXPLICITLY REDUCE BANDED MATRIX TO AN UPPER BANDED MATRIX.
      DO 1 I=J1,J2
      H=X(I+1)-X(I)
      Y0=Y(I+1)-Y(I)
      H2=H*H
      HA=H-A0
      H2A=H-2.D0*A0
      H3A=2.D0*H-3.D0*A0
      H2B=H2*B0
      Q(1,I)=H2/HA
      Q(2,I)=-HA/(H2A*H2)
      Q(3,I)=-H*H2A/H3A
      F(1,I)=(Y0-H*B0)/(H*HA)
      F(2,I)=(H2B-Y0*(2.D0*H-A0))/(H*H2*H2A)
      F(3,I)=-(H2B-3.D0*Y0*HA)/(H*H3A)
      A0=Q(3,I)
 1    B0=F(3,I)
C   TAKE CARE OF LAST TWO ROWS.
      I=J2+1
      H=X(I+1)-X(I)
      Y0=Y(I+1)-Y(I)
      H2=H*H
      HA=H-A0
      H2A=H*HA
      H2B=H2*B0-Y0*(2.D0*H-A0)
      Q(1,I)=H2/HA
      F(1,I)=(Y0-H*B0)/H2A
      HA=X(J2)-X(I+1)
      Y0=-H*HA*(HA+H)
      HA=HA*HA
C   CALCULATE DERIVITIVE AT FAR END.
      Y0=(Y(I+1)*(H2-HA)+Y(I)*HA-Y(J2)*H2)/Y0
      Q(3,I)=(Y0*H2A+H2B)/(H*H2*(H-2.D0*A0))
      Q(2,I)=F(1,I)-Q(1,I)*Q(3,I)
C   SOLVE UPPER BANDED MATRIX BY REVERSE ITERATION.
      DO 2 J=J1,J2
      K=I-1
      Q(1,I)=F(3,K)-Q(3,K)*Q(2,I)
      Q(3,K)=F(2,K)-Q(2,K)*Q(1,I)
      Q(2,K)=F(1,K)-Q(1,K)*Q(3,K)
 2    I=K
      Q(1,I)=B1
C   FILL IN THE LAST POINT WITH A LINEAR EXTRAPOLATION.
 9    J2=J2+2
      DO 14 J=1,3
 14   Q(J,J2)=YY(J)
C   SEE IF THIS DISCONTINUITY IS THE LAST.
 12   IF(J2-I2)6,13,13
C   NO.  GO BACK FOR MORE.
 6    J1=J2+2
      IF(J1-I2)8,8,7
C   THERE IS ONLY ONE POINT LEFT AFTER THE LATEST DISCONTINUITY.
 7    DO 16 J=1,3
 16   Q(J,I2)=YY(J)
C   FINI.
 13   RETURN
      END

c================================================================
      subroutine rsoinc(A,N,IDX)
	implicit real*8 (a-h,o-z)
      DIMENSION A(1),IDX(1) 
c
c	implicit real*8 (a-h,o-z)
c
      IF (N.EQ.1) GO TO 65
      IF (N.LE.0) GO TO 60
      DO 1 I = 1,N
      IDX(I) = I
    1 CONTINUE
      N2 = N/2
      N21 = N2 + 2
      ICT=1 
      I=2 
   11 N1=N21-I
      NN=N
      IK=N1 
   15 C=A(IK) 
      IC=IDX(IK)
  100 JK=2*IK 
      IF (JK.GT.NN) GO TO 140 
      IF (JK.EQ.NN) GO TO 120 
       IF (A(JK+1).LE.A(JK)) GO TO 120
      JK=JK+1 
  120 IF (A(JK).LE. C) GO TO 140
      A(IK)=A(JK) 
      IDX(IK)=IDX(JK) 
      IK=JK 
      GO TO 100 
  140 A(IK)=C 
      IDX(IK)=IC
      GO TO (3,45) ,ICT 
    3 IF (I.GE.N2) GO TO 35 
      I=I+1 
      GO TO 11
   35 ICT=2 
      NP2=N+2 
      I=2 
   37 N1=NP2-I
      NN=N1 
      IK=1
      GO TO 15
  45  CONTINUE
      T = A(1)
      A(1) = A(N1)
      A(N1) = T 
      IT = IDX(1) 
      IDX(1) = IDX(N1)
      IDX(N1) = IT
      IF (I.GE.N) GO TO 55
      I=I+1 
      GO TO 37
   55 RETURN
60	continue
c   60 WRITE(16,500)
  500 FORMAT('ERROR RETURN FROM SORTD1 - N LESS THAN OR EQUAL TO 1')
      STOP
   65 IDX(1)=1
      RETURN
      END


c================================================================
      subroutine extent(NSQ,XLAMIN,XLAMAX,XLOMIN,XLOMAX,
     &	nsqrs,nsqtot,nlatzones,eq_incr_in)
c--finds coordinate range of square number nsq
	IMPLICIT REAL*8 (A-H,O-Z)
	real*4 eq_incr_in
      dimension nsqrs(nlatzones),nsqtot(nlatzones+1)
	eq_incr=eq_incr_in
      LAZONE=2
      DO WHILE (NSQ.GT.NSQTOT(LAZONE))
         LAZONE=LAZONE+1
      ENDDO
      LAZONE=LAZONE-1
      NNSQ=NSQ-NSQTOT(LAZONE)
      XLAMIN=90.-LAZONE*eq_incr
      XLAMAX=XLAMIN+eq_incr
      GRSIZE=360./NSQRS(LAZONE)
      XLOMAX=NNSQ*GRSIZE
      XLOMIN=XLOMAX-GRSIZE
      RETURN
      END

      DOUBLE PRECISION FUNCTION DDABS(X)
      DOUBLE PRECISION X
      IF(X.GE.0.D0) DDABS=X
      IF(X.LT.0.D0) DDABS=-X
      RETURN
      END

c================================================================
c================================================================
