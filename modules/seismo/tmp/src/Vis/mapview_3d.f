
c-----from a 3d mantle model stored in lapo's format, to a set of map-files (one per
c-----model layer) in a format compatible with GMT command psxy
	program points
c	parameter(nlayers=15)
	parameter(nlayersm=60)
c	parameter(nlatzomax=180,iswit=1)
	parameter(nlatzomax=180)
	DIMENSION NSQRS(nlatzomax),ave(nlayersm),nsqtot(nlatzomax+1)
	character name*80,filename3*11,answ*1,answ2*1
	character colorname*30,answl*1
	integer irgb(100,3),rgb(3)
	dimension v(100)
	character*80 string

	print*,"horizontal size of voxels?"
	read*,eq_incr
	print*,"compatibility with reference grid 1=yes 0=no"
	read*,iswit
	refgrid=eq_incr
	nlatzones=180./eq_incr
	print*,'number of layers?'
	read*,nlayers
	if(nlayers.gt.nlayersm)then
	   print*,'out of bounds',nlayersm
	   stop
	endif

	print*,'model file?'
	read*,name

	print*,'palette file?'
	read*,colorname

	print*,'remove average from each layer?'
	read*,answ2
	do i=1,nlayers
	   ave(i)=0.
	enddo

        print*,"top and bottom layers?"
        read*,rtop,rbot
c	layer spacing
	rincr=(rtop-rbot)/nlayers

	print*,'mode (1=GMT layers, 2=values layer, 3=values xyz 4=center of cells)'
	read*,imode
	
        call param(eq_incr,nsqrs,nsqtot,nlatzones,n1layer,iswit,refgrid,nlatzomax)

	if(imode.ne.4)then
	   open(1,file=name,status='old')
	   
	   if(answ2.eq.'y'.or.answ2.eq.'Y')then
	      call layav(nlayers,ave,n1layer,rtop,rbot)
	   endif
	endif

	if(imode.eq.1)then
	   open(unit=54,file=colorname,status='old')
c--read in discrete gmt color palette table 
	   i=1
 1	   read(54,"(a80)",end=2)string
	   if((string(1:1).eq."#").or.(string(1:1).eq."B").or.
     &(string(1:1).eq."F").or.(string(1:1).eq."N"))goto1
	   read(string,*)v(i),(irgb(i,k),k=1,3),v(i+1),(irgb(i+1,k),k=1,3)
	   write(*,*)i,v(i),(irgb(i,k),k=1,3),v(i+1)
	   i=i+1
	   goto 1
 2	   continue 
	endif
	nint=i-1
	if(imode.eq.3)then
	   open(unit=3,file='valxyz.dat')
	   print *,'writing to valxyz.dat'
	else if(imode.eq.4)then
	   open(unit=3,file='centxy.dat')
	   print *,'writing to centxy.dat'
	endif

	if(imode.eq.4)then 
	   nend = 1
	else
	   nend = nlayers
	endif
	do nl=1,nend
	   if(imode.lt.3)then
	      write(filename3,'(a6,i2.2,a2)')'layer_',nl
	      open(unit=3,file=filename3)
	   endif
	   zdepth = rincr*(float(nl)-.5)

	   do ila=1,nlatzones
	      rlati=90.-(eq_incr*(ila-1))
	      do isq=1,nsqrs(ila)

		 rlong=(360./nsqrs(ila))*(isq-1)
		 RINLO=(360./nsqrs(ila))
		 xlomin=RLONG
		 xlomax=RLONG+RINLO
		 xlamax=RLATI
		 xlamin=RLATI-eq_incr
		 if(imode.ne.4)then
c       read value and remove mean
		    read(1,*)k,ccc
		    ccc=(ccc-ave(nl))
		 endif
		 if(imode.eq.1)then
c--     determine color
		    do j=1,nint
		       if((ccc.ge.v(j)).and.(ccc.le.v(j+1)))then
			  do k=1,3
			     a=(irgb(j,k)-irgb(j+1,k))/(v(j)-v(j+1))
			     b=irgb(j,k)-(a*v(j))
			     
c	print*,k,a,b,ccc,rgb(k)
			     
			     
			     rgb(k)=a*ccc+b
c---------------------------------------no interpolation of colours
c       rgb(k)=irgb(j,k)
			  enddo
			  
c	print*,isq,j,v(j),ccc,v(j+1),rgb(1),rgb(2),rgb(3)
c	pause
			  
		       endif
		    enddo
		    if(ccc.ge.v(nint+1))then
		       do k=1,3
			  rgb(k)=irgb(nint+1,k)
		       enddo
		    endif
		    if(ccc.le.v(1))then
		       do k=1,3
			  rgb(k)=irgb(1,k)
		       enddo
		    endif
c       print*,ccc,rgb(1),rgb(2),rgb(3)
 40		    write(3,421)'> -G',rgb(1),'/',rgb(2),'/',rgb(3)
		    write(3,*)xlomin,xlamin
		    write(3,*)xlomin,xlamax
		    write(3,*)xlomax,xlamax
		    write(3,*)xlomax,xlamin
		 else if(imode.eq.2)then
c       output of values only
		    write(3,*)(xlomin+xlomax)/2.,(xlamin+xlamax)/2.,ccc
		 else if(imode.eq.3)then
c       x y z val
		    write(3,*)(xlomin+xlomax)/2.,(xlamin+xlamax)/2.,zdepth,ccc
		 else if(imode.eq.4)then
c       x y locations only
		    write(3,*)(xlomin+xlomax)/2.,(xlamin+xlamax)/2.
		    
		 endif

	      enddo
	   enddo
	   if(imode.lt.3)then
	      close(3)
	   endif
	enddo
	if(imode.ge.3)then
	   close(3)
	endif
	if(imode.ne.4)then
	   close(1)
	endif
421	format(a4,i3.3,a1,i3.3,a1,i3.3)
	end

c----------------------------------------------------------------
	subroutine layav(nlayers,ave,n1layer,rtop,rbot)
	dimension ave(nlayers)
	open(17,file="vprofile.txt")
	open(18,file="dvprofile.txt")
	rincr=(rtop-rbot)/nlayers
	do l=1,nlayers
	   tot=0.
	   tot2=0.
	   do i=1,n1layer
	      read(1,*)ii,value
	      tot=tot+value
	      tot2=tot2+value**2
	      icheck=(n1layer*(l-1))+i
c	      if(ii.ne.icheck)stop
	   enddo

	   std = sqrt ((n1layer * tot2 - tot**2) / ((n1layer*(n1layer-1))));
	   ave(l)=tot/float(n1layer)
c       
	   print*,'a: ',ave(l),(rincr*(float(l)-.5)),rincr*float(l-1),rincr*float(l),std
	   write(17,*)">"
	   write(17,*)rincr*float(l-1),ave(l)
	   write(17,*)rincr*float(l),ave(l)
	   write(18,*)">"
	   write(18,*)rincr*float(l-1),std
	   write(18,*)rincr*float(l),std
	enddo
	close(17)
	close(18)
	print*,'I computed all the ',nlayers,' averages'
	rewind(1)
	return
	end
