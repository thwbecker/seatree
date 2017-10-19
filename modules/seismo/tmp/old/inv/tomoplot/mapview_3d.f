c-----from a 3d mantle model stored in lapo's format, to a set of map-files (one per
c-----model layer) in a format compatible with GMT command psxy
	program points
	parameter(nlayers=15)
c	parameter(eq_incr=10.)
c	parameter(nlatzones=180/eq_incr,refgrid=eq_incr,iswit=0)
	parameter(nlatzomax=180,iswit=0)
c	DIMENSION NSQRS(nlatzones),ave(nlayers),nsqtot(nlatzones+1)
	DIMENSION NSQRS(nlatzomax),ave(nlayers),nsqtot(nlatzomax+1)
	character name*80,filename3*11,answ*1,answ2*1
	character colorname*30,answl*1
	integer irgb(100,3),rgb(3)
	dimension v(100)

	print*,"horizontal size of voxels?"
	read*,eq_incr
	nlatzones=180./eq_incr

	do i=1,nlayers
	   ave(i)=0.
	enddo
	call param(eq_incr,nsqrs,nsqtot,nlatzones,n1layer,iswit,refgrid,
     &	nlatzomax)
cTEST
c	print*,nsqrs
c	print*,n1layer

	print*,'select file to read from'
	read*,name
	print*,'dummy line?(y/n)'
	read*,answl
	ksign=1
	print*,'switch sign?(y/n)'
	read*,answ
	if(answ.eq.'y'.or.answ.eq.'Y')ksign=-1
	open(1,file=name,status='old')
	print*,'what .cpt file?'
	read*,colorname
	open(unit=54,file=colorname,status='old')
	print*,'remove average from each layer?'
	read*,answ2
	if(answ2.eq.'y'.or.answ2.eq.'Y')then
	   call layav(nlayers,ave,answl,n1layer)
	endif
	if(answl.eq.'y'.or.answl.eq.'Y')then
	   read(1,*)
	endif

c--read in discrete gmt color palette table 
	i=1
1	read(54,*,end=2)v(i),(irgb(i,k),k=1,3),v(i+1)
	i=i+1
	goto1
2	nint=i-1

	do 100 nl=1,nlayers
	write(filename3,'(a6,i2.2,a2)')'layer_',nl
	open(unit=3,file=filename3)
	do 500 ila=1,nlatzones
	   rlati=90.-(eq_incr*(ila-1))
	   do 400 isq=1,nsqrs(ila)
	      rlong=(360./nsqrs(ila))*(isq-1)
	      RINLO=(360./nsqrs(ila))
	      xlomin=RLONG
	      xlomax=RLONG+RINLO
	      xlamax=RLATI
	      xlamin=RLATI-eq_incr
	      read(1,*)k,ccc
	      ccc=(ccc-ave(nl))*ksign
c--determine color
	      do j=1,nint
	         if((ccc.ge.v(j)).and.(ccc.lt.v(j+1)))then
	            do k=1,3
                       a=(irgb(j,k)-irgb(j+1,k))/(v(j)-v(j+1))
                       b=irgb(j,k)-(a*v(j))
                       rgb(k)=a*ccc+b
c---------------------------------------no interpolation of colours
c	               rgb(k)=irgb(j,k)
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
40	      write(3,421)'> -G',rgb(1),'/',rgb(2),'/',rgb(3)
	      write(3,*)xlomin,xlamin
	      write(3,*)xlomin,xlamax
	      write(3,*)xlomax,xlamax
	      write(3,*)xlomax,xlamin
400	   continue
500	continue
	close(3)
100	continue
	close(1)
421	format(a4,i3.3,a1,i3.3,a1,i3.3)
	end

c----------------------------------------------------------------
	subroutine layav(nlayers,ave,answl,n1layer)
	dimension ave(nlayers)
	character*1 answl
	open(17,file="vprofile.txt")
	rincr=(6371.-3471.)/nlayers
	if(answl.eq.'y'.or.answl.eq.'Y')then
	   read(1,*)
	endif
	do l=1,nlayers
	   tot=0.
	   do i=1,n1layer
	      read(1,*)ii,value
	      tot=tot+value
	      icheck=(n1layer*(l-1))+i
c	      if(ii.ne.icheck)stop
	   enddo
	   ave(l)=tot/float(n1layer)
	   print*,ave(l),(rincr*(float(l)-.5)),rincr*float(l-1),rincr*float(l)
	   write(17,*)">"
	   write(17,*)rincr*float(l-1),ave(l)
	   write(17,*)rincr*float(l),ave(l)
	enddo
	close(17)
	print*,'I computed all the ',nlayers,' averages'
	rewind(1)
	return
	end
