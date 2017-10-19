	parameter(nlaym=15)
	parameter(n0max=20000)
	parameter(nmax=n0max*nlaym)
	PARAMETER(NONZ=60000000)
	PARAMETER(M=800000)
	parameter(nlatzomax=180)
	parameter(iswit=0)
	dimension indx(nonz),values(nonz),mpoin(0:m)
        integer nsqrs(nlatzomax),nsqtot(nlatzomax+1)
        character*80 filemod
        dimension xmap(nmax)

	print*,"what pixel size?"
	read*,eq_incr
	nlatzones=180/eq_incr
        refgrid=eq_incr
        print*,"what model?"
        read*,filemod

        call param(eq_incr,nsqrs,nsqtot,nlatzones,n1layer,iswit,refgrid,
     &       nlatzomax)
	n_0=n1layer
        nlay=nlaym
	n=n_0*nlay
c	n_1=n_0*nlaym
	if(n_0.gt.n0max)stop "too many pixels"
	if(n.gt.nmax)stop "too many voxels"
c----------------read in model
	print*,"opening",filemod
        open(1,file=filemod,status="old")
        read(1,*)
        do i=1,n
           read(1,*)k,xmap(i)
           if(k.ne.i)stop "wrong file format"
        enddo
        close(1)
c----------------find roughness damping matrix
        inp=0
        nelm=0
        nelrhs=0
        call gradamp4(1.,1.,nelm,nelrhs,
     &       values,indx,mpoin,
     &       m,n,NONZ,nlaym,inp,nsqrs,nsqtot,nlatzones,n_0,eq_incr,
     &       nlatzomax)
c----------------dot product with model vector
        roughness=0.
        do i=1,nelrhs
           rowprod=0.
           do j=mpoin(i-1)+1,mpoin(i)
              rowprod=rowprod+(values(j)*xmap(indx(j)))
           enddo
           roughness=roughness+(rowprod*rowprod)
        enddo
        do k=1,80
           if(filemod(k:k).eq." ")goto1
        enddo
 1      kname=k-1
        open(424,file="roughness.txt",access="append")
        write(424,*)filemod(1:kname),roughness
        close(424)
        print*,filemod(1:kname),roughness
        end
c------------------------------------------------------------------------
	subroutine gradamp4(weight,weightv,nnn,nelp,xxx,indx,mpoin,m,n,
     &nonz,nlay,koffs,nsqrs,nsqtot,nlatzones,n_0,eq_incr,
     &nlatzomax)
c--defines the matrix corresponding to gradient damping
c--and appends it to the A matrix, directly in memory.
c	dimension nsqrs(nlatzones),nsqtot(nlatzones+1)
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
	dimension indx(nonz),xxx(nonz),mpoin(0:m) !,t(m)
	parameter(zer=0.)

	DO 1000 L=1,NLAY
	more=n_0*(L-1)+KOFFS
        vw=2.
	if((l.eq.1).or.(l.eq.nlay))vw=1.

	do i1=1,nlatzones
	   rloin_2=180./nsqrs(i1)
	   ifirst=nsqtot(i1)+1
	   ilast=nsqtot(i1+1)
	   do i2=ifirst,ilast
c--first, find ileft and iright
	      ileft=i2-1
	      if(i2.eq.ifirst)ileft=ilast
	      iright=i2+1
	      if(i2.eq.ilast)iright=ifirst
c--then, iup and idw...
	      IF((1.LE.I2).AND.(I2.LE.nsqrs(1)))THEN
c	         call coordsurf(i2,rila,rilo)
	         call coordsuper(nlatzomax,i2,rila,rilo,nsqrs,nlatzones,eq_incr)
	         dwla=rila-eq_incr
	         xlol=rilo-rloin_2
	         xlor=rilo+rloin_2
		 idwl=isqre(dwla,xlol,nsqrs,nsqtot,nlatzones,n_0,eq_incr)
	         idwr=isqre(dwla,xlor,nsqrs,nsqtot,nlatzones,n_0,eq_incr)
c--increment the matrix
	         xxx(nnn+1)=(3.+vw)*weight
	         indx(nnn+1)=i2+MORE
	         xxx(nnn+2)=-weight
	         indx(nnn+2)=iright+MORE
	         xxx(nnn+3)=-weight
	         indx(nnn+3)=ileft+MORE
	         nnn=nnn+3
c*****************
	         do idw=idwl,idwr
	         xxx(nnn+1)=-weight/(idwr-idwl+1)
	         indx(nnn+1)=idw+MORE
	         nnn=nnn+1
	         enddo
c*****************
C====================================================
	      elseIF((nsqrs(1)+1.LE.I2).AND.(I2.LE.n_0-nsqrs(nlatzones)))THEN
                 call coordsuper(nlatzomax,i2,rila,rilo,nsqrs,nlatzones,eq_incr)
	         dwla=rila-eq_incr
	         upla=rila+eq_incr
	         xlol=rilo-rloin_2
	         xlor=rilo+rloin_2
	         idwl=isqre(dwla,xlol,nsqrs,nsqtot,nlatzones,n_0,eq_incr)
	         idwr=isqre(dwla,xlor,nsqrs,nsqtot,nlatzones,n_0,eq_incr)
	         iupl=isqre(upla,xlol,nsqrs,nsqtot,nlatzones,n_0,eq_incr)
	         iupr=isqre(upla,xlor,nsqrs,nsqtot,nlatzones,n_0,eq_incr)
c--increment the matrix
	         xxx(nnn+1)=(4.+vw)*weight
	         indx(nnn+1)=i2+MORE
	         xxx(nnn+2)=-weight
	         indx(nnn+2)=iright+MORE
	         xxx(nnn+3)=-weight
	         indx(nnn+3)=ileft+MORE
	         nnn=nnn+3
c*****************
	         do idw=idwl,idwr
	         xxx(nnn+1)=-weight/(idwr-idwl+1)
	         indx(nnn+1)=idw+MORE
	         nnn=nnn+1
	         enddo
c*****************
	         do iup=iupl,iupr
	         xxx(nnn+1)=-weight/(iupr-iupl+1)
	         indx(nnn+1)=iup+MORE
	         nnn=nnn+1
	         enddo
c*****************
C======================================================
	      elseIF((n_0-nsqrs(nlatzones)+1.LE.I2).AND.(I2.LE.n_0))THEN
                 call coordsuper(nlatzomax,i2,rila,rilo,nsqrs,nlatzones,eq_incr)
	         upla=rila+eq_incr
	         xlol=rilo-rloin_2
	         xlor=rilo+rloin_2
	         iupl=isqre(upla,xlol,nsqrs,nsqtot,nlatzones,n_0,eq_incr)
	         iupr=isqre(upla,xlor,nsqrs,nsqtot,nlatzones,n_0,eq_incr)
c--increment the matrix
	         xxx(nnn+1)=(3.+vw)*weight
	         indx(nnn+1)=i2+MORE
	         xxx(nnn+2)=-weight
	         indx(nnn+2)=iright+MORE
	         xxx(nnn+3)=-weight
	         indx(nnn+3)=ileft+MORE
	         nnn=nnn+3
c*****************
	         do iup=iupl,iupr
	         xxx(nnn+1)=-weight/(iupr-iupl+1)
	         indx(nnn+1)=iup+MORE
	         nnn=nnn+1
	         enddo
c*****************
	      ENDIF
c--vertical damping
	   if(l.eq.1)then
	      xxx(nnn+1)=-weightV
	      indx(nnn+1)=i2+n_0+more
	      nnn=nnn+1
	   elseif((l.gt.1).and.(l.lt.nlay))then
	      xxx(nnn+1)=-weightV
	      indx(nnn+1)=i2-n_0+more
	      xxx(nnn+2)=-weightV
	      indx(nnn+2)=i2+n_0+more
	      nnn=nnn+2
	   elseif(l.eq.nlay)then
	      xxx(nnn+1)=-weightV
	      indx(nnn+1)=i2-n_0+more
	      nnn=nnn+1
	   endif

	    nelp=nelp+1
	    mpoin(nelp)=nnn
c	    t(nelp)=zer
	   enddo
	enddo


1000	CONTINUE

	END
c------------------------------------------------------------------------
	function isqre(lat,lon,nsqrs,nsqtot,nlatzones,n,eq_incr)
c----finds the index of the square where (lat,lon) is
	real*4 lat,lon,loc_incr
	dimension nsqrs(nlatzones),nsqtot(nlatzones+1)
	lazone=(90.-lat)/eq_incr+1
	if((90.-lat).gt.180.)lazone=nlatzones
	if((90.-lat).gt.181.)stop "problems in function isqre"
	if(lazone.gt.nlatzones)then
	   print*,"problems in function isqre, latitude",lazone,lat
	   stop
	endif
	if(lon.lt.0.)lon=360.+lon
	if(lon.eq.360.)lon=0.
	loc_incr=360./float(nsqrs(lazone))
	isqre=(lon/loc_incr)+1
	isqre=isqre+nsqtot(lazone)
	if(isqre.gt.n)then
	   print*,"problems in function isqre, longitude"
	   stop
	endif
	RETURN
	END
c------------------------------------------------------------------------
        subroutine param(eq_incr,nsqrs,nsqtot,nlatzones,numto,iswit,
     &refgrid,nlatzomax)
c---find vectors nsqrs and nsqtot that define a block parameterization
c---eq_incr,iswit,refgrid are input, the rest is output
c        dimension nsqrs(nlatzones),nsqtot(nlatzones+1)
        dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
	parameter(pi=3.1415926536)
	numto=0
	colat=-eq_incr/2.
	do k=1,nlatzones
	   colat=colat+eq_incr
	   theta=(colat/180.)*pi
c--for this latitudinal zone, compute number of blocks (nsqrs)
	   deltalon=eq_incr/(sin(theta))
	   nsqrs(k)=(360./deltalon)+1
           if(mod(nsqrs(k),2).ne.0)nsqrs(k)=nsqrs(k)-1
c--if requested, correct nsqrs(k) so the grid is compatible to reference grid
	   if(iswit.eq.1)then
              if(360./nsqrs(k).ge.refgrid)then
 100             if(mod(360./nsqrs(k),refgrid).ne.0)then
                    nsqrs(k)=nsqrs(k)+1
                    goto100
                 else
                 endif
              elseif(360./nsqrs(k).lt.refgrid)then
 101             if(mod(refgrid,360./nsqrs(k)).ne.0)then
                    nsqrs(k)=nsqrs(k)-1
                    goto101
                 else
                 endif
              endif
	   endif
           if(MOD(NSQRS(K),2).ne.0)stop "nsqrs has to be even"
	   nsqtot(k)=numto
	   numto=numto+nsqrs(k)
	enddo
	nsqtot(nlatzones+1)=numto
	print*,nsqrs
	print*,"total number of blocks:",nsqtot(nlatzones+1)
	return
	end
c------------------------------------------------------------------------
	subroutine coordsuper(nlatzomax,nbloc,blocla,bloclo,nsqrs,nlatzones,eq_incr)
c--given a cell index on the Earth's surface, finds longitude and latitude
c--(NOT colatitude) of its center.
c	DIMENSION NSQRS(nlatzones)
	DIMENSION NSQRS(nlatzomax)

	ntot=0
c--loop(s) over all the blocks
	do 500 ila=1,nlatzones
c--increment latitude
	   rlati=90.-(eq_incr*(ila-1))
c--calculate increment in longitude for this band
	   RINLO=(360./nsqrs(ila))
	   do 400 isq=1,nsqrs(ila)
	      rlong=(360./nsqrs(ila))*(isq-1)
	      ntot=ntot+1
	      if(ntot.eq.nbloc)then
	         bloclo=rlong+(rinlo/2.)
	         blocla=rlati-(eq_incr/2.)
	         goto 600
	      endif
400	   continue
500	continue
600	return
	end
c------------------------------------------------------------------------
