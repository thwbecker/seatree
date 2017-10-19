        include 'common_para.h'
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
	dimension nsqrsh(nlatzhmax),nsqtoth(nlatzhmax+1)
	dimension iresolu(500000),ire2(500000)
	dimension inew(500000),iold(500000),inewh(500000),ioldh(500000)
	character filename*80,dummy*80,outpo*80,colorname*80
	integer irgb(100,3),rgb(3)
	dimension v(100)

c--open input and output files
	print*,'file to be read?'
	read*,filename
	open(unit=1,file=filename,status='old')
	print*,'ouptput file?'
	read*,outpo
	open(unit=3,file=outpo)
c	read(1,*)dummy
	print*,'what .cpt file?'
	read*,colorname
	open(unit=54,file=colorname,status='old')
	print*,"size of coarse pixels?"
	read*,eq_incr
	refgrid=eq_incr*1.
	nlatzones=180./eq_incr
	if(nlatzones.gt.nlatzomax)stop "coarse pixels too small"
	print*,"ratio to size of fine pixels?"
	read*,ifa
	eq_hi=eq_incr/ifa
	nlatzohi=180./eq_hi
	if(nlatzohi.gt.nlatzhmax)stop "fine pixels too fine"
 18     print*,"xyz output: coarse (1) or fine (2) grid?"
        read*,ixyz
        if(ixyz.ne.1.and.ixyz.ne.2)goto18
        do k=1,80
           if(outpo(k:k).eq." ")goto19
        enddo
 19     open(18,file=outpo(1:k-1)//".xyz")

c--read in the color palette table
	i=1
1	read(54,*,end=2)v(i),(irgb(i,k),k=1,3),v(i+1),(irgb(i+1,k),k=1,3)
	i=i+1
	goto1
2	nint=i-1

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
	   do j=1,nint
	      if((ccc.ge.v(j)).and.(ccc.lt.v(j+1)))then
	         do k=1,3
	            a=(irgb(j,k)-irgb(j+1,k))/(v(j)-v(j+1))
	            b=irgb(j,k)-(a*v(j))
	            rgb(k)=a*ccc+b
	         enddo
	         goto40	
              elseif(ccc.ge.v(nint+1))then
	         do k=1,3
	            rgb(k)=irgb(nint+1,k)
	         enddo
	         goto40	
	      elseif(ccc.lt.v(1))then
	         do k=1,3
	            rgb(k)=irgb(nint+1,k)
	         enddo
	         goto40	
	      endif
	   enddo
40	   write(3,421)'> -G',rgb(1),'/',rgb(2),'/',rgb(3)
ctest
c     if(xlomin.eq.0.)xlomin=0.001
           if(xlomax.eq.360.)xlomax=359.99
	   write(3,*)xlomin,xlamin
	   write(3,*)xlomin,xlamax
	   write(3,*)xlomax,xlamax
	   write(3,*)xlomax,xlamin
c----------------------also write map in xyz format on another file
           celo=(xlomin+xlomax)/2.
           cela=(xlamin+xlamax)/2.
           if(ixyz.eq.1)write(18,*)celo,cela,ccc
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
	   do j=1,nint
	      if((ccc.ge.v(j)).and.(ccc.lt.v(j+1)))then
	         do k=1,3
	            a=(irgb(j,k)-irgb(j+1,k))/(v(j)-v(j+1))
	            b=irgb(j,k)-(a*v(j))
	            rgb(k)=a*ccc+b
	         enddo
	         goto41	
              elseif(ccc.ge.v(nint+1))then
	         do k=1,3
	            rgb(k)=irgb(nint+1,k)
	         enddo
 	         goto41	
             elseif(ccc.lt.v(1))then
	         do k=1,3
	            rgb(k)=irgb(nint+1,k)
	         enddo
	         goto41	
	      endif
	   enddo
41	   continue
	   write(3,421)'> -G',rgb(1),'/',rgb(2),'/',rgb(3)
	if(xlomax.eq.360.)xlomax=359.99
	   write(3,*)xlomin,xlamin
	   write(3,*)xlomin,xlamax
	   write(3,*)xlomax,xlamax
	   write(3,*)xlomax,xlamin
c----------------------also write map in xyz format on another file
           celo=(xlomin+xlomax)/2.
           cela=(xlamin+xlamax)/2.
           if(ixyz.eq.2)write(18,*)celo,cela,ccc
c--print out a table of the block indexes/locations
c	   write(4,*)ksav,xlomin,xlamin,xlomax,xlamax
	enddo
	write(3,'(a)')'>'

421	format(a4,i3.3,a1,i3.3,a1,i3.3)

	close(1)
	close(3)
	close(54)
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
