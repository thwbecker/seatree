        include 'common_para.h'
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
	dimension nsqrsh(nlatzhmax),nsqtoth(nlatzhmax+1)
	dimension iresolu(500000),ire2(500000)
	dimension inew(500000),iold(500000),inewh(500000),ioldh(500000)
	character filename*200,dummy*200,outpo*200,colorname*200
	integer iclo(100,3),ichi(100,3),rgb(3)
        integer use_cpt
	dimension zlo(100),zhi(100)

c--open input and output files
	print*,'solution file to be read?'
	read*,filename
	open(unit=1,file=filename,status='old')
	print*,'output file?'
	read*,outpo
	open(unit=3,file=outpo)
c	read(1,*)dummy
	print*,'what .cpt file?'
	read*,colorname
        IF (colorname .EQ. 'zout') THEN
           use_cpt = 0
        else
           open(unit=54,file=colorname,status='old')
           use_cpt = 1
        endif
	print*,"size of coarse pixels?"
	read*,eq_incr
	refgrid=eq_incr*1.
	nlatzones=180./eq_incr
	if(nlatzones.gt.nlatzomax)stop "coarse pixels too small"
	print*,"ratio to size of fine pixels?"
	read*,ifa
	if(ifa.ne.1)then
	   print *,'error: this version does not support local refinement'
	   print *,'recompile with proper boundary'
           print *,'settings for high res region'
	   stop
	endif

	eq_hi=eq_incr/ifa
	nlatzohi=180./eq_hi
	if(nlatzohi.gt.nlatzhmax)stop "fine pixels too fine"
        if(use_cpt.eq.1)then
c--   read in the color palette table: one slice per record,
c--   lower z, lower r g b, upper z, upper r g b
           i=1
 1         read(54,*,end=2)zlo(i),(iclo(i,k),k=1,3),
     &          zhi(i),(ichi(i,k),k=1,3)
           i=i+1
           if(i.gt.100)stop 'blk2gmt: more than 99 slices in cpt file'
           goto 1
 2         nint=i-1
        else
           nint=0
        endif
        

c--determine nsqrs,nsqrsh
	numto=0
	numhi=0
	colat=-eq_incr/2.
	do k=1,nlatzones
	   colat=colat+eq_incr
	   theta=(colat/180.)*pi
	   deltalon=eq_incr/(sin(theta))
	   nsqrs(k)=(360./deltalon)+1
c----needs to be even
	   if(mod(nsqrs(k),2).ne.0)nsqrs(k)=nsqrs(k)-1
c-------------------------------new
c--if requested, correct nsqrs(k) so the grid is compatible to reference grid
	   if(iswit.eq.1)then
	    if(360./nsqrs(k).ge.refgrid)then
100	     if(mod(360./nsqrs(k),refgrid).ne.0)then
	      nsqrs(k)=nsqrs(k)+1
	      goto100
	     else
	     endif
	    elseif(360./nsqrs(k).lt.refgrid)then
101	     if(mod(refgrid,360./nsqrs(k)).ne.0)then
	      nsqrs(k)=nsqrs(k)+1
	      goto101
	     else
	     endif
	    endif
	   endif
c----------------------------------
	   do j=1,ifa
	      kfine=((k-1)*ifa)+j
	      nsqrsh(kfine)=nsqrs(k)*ifa
	      nsqtoth(kfine)=numhi
	      numhi=numhi+nsqrsh(kfine)
	   enddo
	   nsqtot(k)=numto
	   numto=numto+nsqrs(k)
	enddo
	nsqtot(nlatzones+1)=numto
	nsqtoth(nlatzohi+1)=numhi
	print*,'numto=',numto,'  numhi=',numhi

c--determine iresolu,ire2,kireso,kire2
	kireso=0
	kire2=0
	print*,westbo,eastbo,southbo,rthnobo,eq_incr
	do parall=rthnobo,southbo,-eq_incr
	   do rmerid=westbo,eastbo,eq_incr
	      kireso=kireso+1
	      ilat=parall*100.+0.5
	      ilon=rmerid*100.+0.5
	      iresolu(kireso)=
     &isqre(parall,rmerid,nsqrs,nsqtot,nlatzones,numto,eq_incr)
c     &	superisqre(ilat,ilon,nsqrs,nsqtot,nlatzones,numto,eq_incr)
	      icoarse=iresolu(kireso)
	      call rang(icoarse,xlamin,xlamax,
     &	xlomin,xlomax,nsqrs,nsqtot,nlatzones,n,eq_incr)
	      do ifila=1,ifa
	      do ifilo=1,ifa
	         kire2=kire2+1
	         xlafi=xlamin+((xlamax-xlamin)/ifa)*(ifila-0.5)
	         xlofi=xlomin+((xlomax-xlomin)/ifa)*(ifilo-0.5)
	         ilat=xlafi*100.+0.5
	         ilon=xlofi*100.+0.5
c	         ire2(kire2)=
c     &	superisqre(ilat,ilon,nsqrsh,nsqtoth,nlatzohi,numhi,eq_hi)
	         ire2(kire2)=
     &	isqre(xlafi,xlofi,nsqrsh,nsqtoth,nlatzohi,numhi,eq_hi)
	      enddo
	      enddo
	   enddo
	enddo
	print*,kireso,' nonzero elements in array iresolu'
	print*,kire2,' nonzero elements in array ire2'

c--count blocks used in parameterization
	icoar=0
	do icoblo=1,numto
	   do iche=1,kireso
	      if(icoblo.eq.iresolu(iche))then
	         icoar=icoar+1
	         goto37
	      endif
	   enddo
37	continue
	enddo
	ifine=0
	do icoblo=1,numhi
	   do iche=1,kire2
	      if(icoblo.eq.ire2(iche))then
	         ifine=ifine+1
	         goto38
	      endif
	   enddo
38	continue
	enddo
	nprime=numto-icoar+ifine
	print*,"nprime=",nprime

c--determine inew,inewh,iold,ioldh
	call correspc(iresolu,kireso,numto,0,inew,iold,ico)
	print*,'compare:',numto-icoar,inew(numto)
	iof=numto-icoar
	call corresph(ire2,kire2,numhi,iof,inewh,ioldh,ico)

c--read the coarse part of the model and augment the map file
	do i=1,numto-icoar
	   read(1,*)k,ccc
ctest
c	print*,k

	   indblo=iold(k)
	   ksav=k
	   call rang(indblo,xlamin,xlamax,
     &	xlomin,xlomax,nsqrs,nsqtot,nlatzones,numto,eq_incr)
           if(use_cpt.eq.1)then 
              call cptcol(ccc,nint,zlo,zhi,iclo,ichi,rgb)
 40           write(3,421)'> -G',rgb(1),'/',rgb(2),'/',rgb(3)
           else
              write(3,*)'> -Z',ccc
           endif
ctest
c	if(xlomin.eq.0.)xlomin=0.001
	if(xlomax.eq.360.)xlomax=359.99
	   write(3,*)xlomin,xlamin
	   write(3,*)xlomin,xlamax
	   write(3,*)xlomax,xlamax
	   write(3,*)xlomax,xlamin
c--print out a table of the block indexes/locations
c	   write(4,*)ksav,xlomin,xlamin,xlomax,xlamax
	enddo

ctest
	print*,"blk2gmt okay"

c--read the fine part of the model and augment the map file
	do i=1,ifine
	   read(1,*)k,ccc
	   indblo=ioldh(k)
	   ksav=k
	   call rang(indblo,xlamin,xlamax,
     &	xlomin,xlomax,nsqrsh,nsqtoth,nlatzohi,numhi,eq_hi)
           if(use_cpt.eq.1)then
              call cptcol(ccc,nint,zlo,zhi,iclo,ichi,rgb)
              write(3,421)'> -G',rgb(1),'/',rgb(2),'/',rgb(3)
           else
              write(3,*)'> -Z',ccc
           endif
	if(xlomax.eq.360.)xlomax=359.99
	   write(3,*)xlomin,xlamin
	   write(3,*)xlomin,xlamax
	   write(3,*)xlomax,xlamax
	   write(3,*)xlomax,xlamin
c--print out a table of the block indexes/locations
c	   write(4,*)ksav,xlomin,xlamin,xlomax,xlamax
	enddo
	write(3,'(a)')'>'

421	format(a4,i3.3,a1,i3.3,a1,i3.3)

	close(1)
	close(3)
        if(use_cpt.eq.1)then
           close(54)
        endif
	end

c************************************************************
	subroutine correspc(iresolu,kireso,n,ioffset,inew,iold,ico)
	dimension iresolu(500000),inew(500000),iold(500000)
	ico=0
	do i=1,n
	   do k=1,kireso
	   if(i.eq.iresolu(k))then
	      ico=ico+1
	      inew(i)=-1
	      goto42
	   endif
	   enddo
	   inew(i)=(i+ioffset)-ico
42	   iold(inew(i))=i
	enddo
	return
	end
c************************************************************
	subroutine corresph(iresolu,kireso,n,ioffset,inew,iold,ico)
	dimension iresolu(500000),inew(500000),iold(500000)
	ico=0
	do i=1,n
	   do k=1,kireso
	   if(i.eq.iresolu(k))then
	      inew(i)=(i+ioffset)-ico
	      goto42
	   endif
	   enddo
	   ico=ico+1
	   inew(i)=-1
42	   iold(inew(i))=i
	enddo
	return
	end
c************************************************************
	function isqre(xlat,xlon,nsqrs,nsqtot,nlatzones,n,eq_incr)
c--finds the number of the square where (xlat,xlon) is

c	dimension nsqrs(nlatzones),nsqtot(nlatzones+1)
	parameter(nlatzomax=180,nlatzhmax=720)
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
	dimension nsqrsh(nlatzhmax),nsqtoth(nlatzhmax+1)
       	lazone=(90.-xlat)/eq_incr+1
	if(lazone.gt.nlatzones)lazone=nlatzones
c	llon=lon
c	if(llon.lt.0)llon=36000+llon
	isqre=(xlon/360.)*nsqrs(lazone)+1
	isqre=isqre+nsqtot(lazone)
c	if(isqre.gt.n)isqre=n
	return
	end
c************************************************************
	subroutine rang(nsq,xlamin,xlamax,xlomin,xlomax,
     &	nsqrs,nsqtot,nlatzones,n,eq_incr)
c
c finds the coordinate range of square number 'nsq'
c
c	dimension nsqrs(nlatzones),nsqtot(nlatzones+1)
	parameter(nlatzomax=180,nlatzhmax=720)
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)

	lazone=2
	do while (nsq.gt.nsqtot(lazone))
	   lazone=lazone+1
	enddo
	lazone=lazone-1
	nnsq=nsq-nsqtot(lazone)
	xlamin=90.-lazone*eq_incr
	xlamax=xlamin+eq_incr
	grsize=360./nsqrs(lazone)
	xlomax=nnsq*grsize
	xlomin=xlomax-grsize
	return
	end
c************************************************************

c
c     color for value ccc from a color palette table with ncpt slices
c     (lower bound zlo, upper bound zhi, colors iclo and ichi at those
c     bounds). Within a slice, the color is interpolated linearly
c     between the lower and upper color, which returns the slice color
c     for discrete tables and the interpolated color for continuous
c     ones. Values below or above the table get the lowest or highest
c     color.
c
      subroutine cptcol(ccc,ncpt,zlo,zhi,iclo,ichi,rgb)
      real ccc,zlo(100),zhi(100)
      integer ncpt,iclo(100,3),ichi(100,3),rgb(3)
      integer j,k,jsel
      real frac
      if(ccc.lt.zlo(1))then
         do k=1,3
            rgb(k)=iclo(1,k)
         enddo
         return
      endif
      if(ccc.ge.zhi(ncpt))then
         do k=1,3
            rgb(k)=ichi(ncpt,k)
         enddo
         return
      endif
c     find the slice; if the value falls into a gap between slices,
c     the last slice below it is used
      jsel=1
      do j=1,ncpt
         if(ccc.ge.zlo(j))then
            jsel=j
         endif
      enddo
      if(zhi(jsel).gt.zlo(jsel))then
         frac=(ccc-zlo(jsel))/(zhi(jsel)-zlo(jsel))
      else
         frac=0.
      endif
      if(frac.gt.1.)then
         frac=1.
      endif
      do k=1,3
         rgb(k)=nint(iclo(jsel,k)+frac*(ichi(jsel,k)-iclo(jsel,k)))
      enddo
      return
      end
