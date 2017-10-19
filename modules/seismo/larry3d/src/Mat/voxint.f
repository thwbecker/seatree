	implicit real*8 (a-h,o-z)
	parameter(nint=4) !how many interpolations of ray path
	parameter(rearth=6371.d0,rcmb=3480.d0) ! careful, this is hardwired
	parameter(nlatzomax=180)
	parameter(ndim=10000)
	parameter(toler=0.00005)! numerical error in pixel coordinates
        character*200 dtfile,qkfile,stfile,modfile
	character chdatin*6 !test
        character*256 raynam,filnam
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
        dimension rbnd(0:60)
        real*4 rescale(60),elat,elon,slat,slon,delta,azep,azst,datum,depth
        common/lapo/raddown(ndim),radup(ndim),kdown,kup,time(ndim),
     &dedown(ndim),deup(ndim),raydelr(ndim),z(ndim),iswtc
	dimension x(ndim),y(ndim),raydel(ndim)
	dimension raydtmp(ndim),ztmp(ndim),idx(ndim)
	dimension xi(ndim*2),yi(ndim*2),delint(ndim*2),zi(ndim*2)
	dimension radmod(ndim),vstart(ndim),work1(3,ndim),work2(3,ndim)
	pi=2.*asin(1.)
	radian=180./pi
	filnam="aniprmc" !reference model, now hardwired

	print*,"size of pixels?"
	read*,eq_incr
	nlatzones=180./eq_incr
	if(nlatzones.gt.nlatzomax)stop "coarse pixels too small"
	print*,"how many layers?"
	read*,nlay
        print*,"what data?"
        read*,raynam
	print*,"minimum acceptable epicentral distance (degrees)?"
	read*,deltamin
3       print*,"polarization? (0 for SV, 1 for SH)"
        read*,ipolar
        if(ipolar.ne.0.and.ipolar.ne.1)goto3
        print*,"file with data?"
        read*,dtfile
        print*,"file with sources?"
        read*,qkfile
        print*,"file with stations?"
        read*,stfile
        open(921,file=dtfile,access='direct',recl=4,status='old',form='unformatted')
        open(422,file=qkfile,access='direct',recl=12,status='old',form='unformatted')
        open(522,file=stfile,access='direct',recl=8,status='old',form='unformatted')
	print*,"file with starting model?"
	read*,modfile
	open(1,file=modfile,status="old")
	i=1
4	read(1,*,end=5)radmod(i),vstart(i)
cTEST
c	write(79,*)radmod(i),vstart(i)
	i=i+1
	goto4
5	nstartmod=i-1
	if(radmod(2).ge.radmod(1))stop "sort starting model file"

c--define nsqrs,nsqtot and layers
	call param(eq_incr,nsqrs,nsqtot,nlatzones,n1layer,0,eq_incr,nlatzomax)
	nvx=n1layer*nlay
	print*,nvx," voxels total" !test
	rbnd(0)=rearth
	yncr=(6371.-rcmb)/dble(nlay)
	do iyn=1,nlay
	   rbnd(iyn)=6371.-(yncr*dble(iyn))
	enddo

c--open the files where the matrix will be stored
        open(111,file='a.vx_mat',access='direct',recl=4,form='unformatted')
        open(112,file='a.vx_ind',access='direct',recl=4,form='unformatted')
        open(113,file='a.vx_pnt')
        open(109,file='d.vx_vec')
        open(50,file="used_sources_stations.txt")

	j=1
	irec=0
c----------read datum
10      read(921,rec=j,err=20)datum
        read(422,rec=j)elon,elat,depth
        read(522,rec=j)slon,slat

cTEST
	if(elat.eq.slat.and.(abs(elon-slon).eq.180.))slat=slat-toler
c	print*,elon,elat,depth,slon,slat

	if(mod(j,1000).eq.0)print*,j," couples processed"
        if(elon.lt.0.)elon=elon+360.
        if(slon.lt.0.)slon=slon+360.
c-----------find epicentral distance
        call delazs(elat,elon,slat,slon,delta,azep,azst)
c---------------------------------exclude short epicentral distances
	if(delta.lt.deltamin)then
	j=j+1
	goto10
	endif

        ddel=dble(delta)
        ddep=dble(depth)
	nto6=0
	nto7=0
c-----------find ray path
	iflapkp=0 ! temporary fix
	call find_path(ddel,raynam,ddep,1,0,filnam,nto6,nto7,iflapkp,ipolar)
	if(nto6.ne.0.or.nto7.ne.0)then
	   print*,j,"failure of find_path"
	   j=j+1
	   goto10
	endif
	eplon=dble(elon)
	eplat=dble(elat)
	stlon=dble(slon)
	stlat=dble(slat)

	do i=1,kup	
	   raydel(i)=deup(i)*radian
	   if(raydel(i).lt.0.)then
	      print*,j,"failure of find_path"
	      j=j+1
	      goto10
	   endif
	   z(i)=radup(i)
c	   write(87,*)deup(i)*radian,radup(i)
	enddo
	do i=1,kdown
	   raydel(kup+i)=dedown(i)*radian
	   if(raydel(kup+i).lt.0.)then
	      print*,j,"failure of find_path"
	      j=j+1
	      goto10
	   endif
	   z(kup+i)=raddown(i)
c	   write(87,*)dedown(i)*radian,raddown(i)
	enddo
	npt=kup+kdown

c--sort ray path points in order of increasing delta
	call rsoinc(raydel,npt,idx)
	do i=1,npt
	   ztmp(i)=z(idx(i))
	enddo
	do i=1,npt
	   z(i)=ztmp(i)
	enddo	

c--find_path routine occasionally repeats points
	k=1
	raydtmp(1)=raydel(1)
	ztmp(1)=z(1)
	do i=2,npt
c	   if(raydel(i).ne.raydel(i-1))then
	   if(abs(raydel(i)-raydel(i-1)).gt.toler)then
	      k=k+1
	      raydtmp(k)=raydel(i)
	      ztmp(k)=z(i)
	   endif
	enddo
	npt=k
	do i=1,npt
	   raydel(i)=raydtmp(i)
	   z(i)=ztmp(i)
	enddo

	call gceq(eplon,eplat,stlon,stlat,raydel,npt,x,y,ndim)

	if(raydel(1).ne.0.)then
	   x(npt+2)=stlon
	   y(npt+2)=stlat
	   z(npt+2)=rearth
	   raydel(npt+2)=ddel
	   do i=npt+1,2,-1
	      y(i)=y(i-1)
	      x(i)=x(i-1)
	      z(i)=z(i-1)
	      raydel(i)=raydel(i-1)
	   enddo
	   x(1)=eplon
	   y(1)=eplat
	   z(1)=rearth-ddep
	   raydel(1)=0.d0
	   npt=npt+2
	else
	   x(npt+1)=stlon
	   y(npt+1)=stlat
	   z(npt+1)=rearth
	   raydel(npt+1)=ddel
	   x(1)=eplon
	   y(1)=eplat
	   z(1)=rearth-ddep
	   raydel(1)=0.d0
	   npt=npt+1
	endif
	
cTEST
	write(chdatin,"(i6.6)")j
c	open(84,file="test/A"//chdatin//".txt")
c	open(99,file="int"//chdatin//".txt")
c	open(98,file="raypath"//chdatin//".txt")
c	do i=1,npt
c	   write(98,"(4(f12.5,2x))")x(i),y(i),z(i),raydel(i)
c	enddo
c	close(98)

	do k=1,nint
	   call linint(y,yi,raydel,delint,ndim,npt)
	   call linint_x(x,xi,raydel,delint,ndim,npt)
	   call linint(z,zi,raydel,delint,ndim,npt)
	   do i=1,npt*2-1
	      if(xi(i).lt.0.)xi(i)=360.+xi(i)
	      if(xi(i).gt.360.)xi(i)=xi(i)-360.
	      x(i)=xi(i)
	      y(i)=yi(i)
	      z(i)=zi(i)
	      raydel(i)=delint(i)
	   enddo
	   npt=2*npt-1
	enddo
	
c	open(98,file="raypath_int"//chdatin//".txt")
c	do i=1,npt
c	   write(98,"(4(f12.5,2x))")x(i),y(i),z(i),raydel(i)
c	enddo
c	close(98)

c---------------------project ray path onto grid
	do ila=1,nlay
	   if((z(1).le.rbnd(ila-1)).and.(z(1).ge.rbnd(ila)))then
	   iv0=ila
	   endif
        enddo
	ih0=isqre(y(1),x(1),nsqrs,nsqtot,nlatzones,n1layer,eq_incr)
	ind0=(iv0-1)*n1layer+ih0
	call span(ih0,ymi0,yma0,xmi0,xma0,nsqrs,nsqtot,nlatzones,eq_incr)
	x0=x(1)
	y0=y(1)
	z0=z(1)
	d0=0.

	do i=1,npt !loop over all points in ray path

cTEST
c	write(*,"(i6,1x,3(f12.6,2x),i6)")i,x0,y0,z0,ind0

c	if(xma0.eq.360.)xma0=0.
c	if(xmi0.eq.0.)xmi0=360.
        if(abs(xma0-360.).lt.toler)xma0=0.
        if(abs(xmi0).lt.toler)xmi0=360.

c--determine index of voxel for i-th point on ray path
	   do ila=1,nlay
              if((z(i).le.rbnd(ila-1)).and.(z(i).ge.rbnd(ila)))then
                 iv=ila
              endif
           enddo
	   if(z(i).lt.rbnd(nlay))iv=nlay+1
	   ih=isqre(y(i),x(i),nsqrs,nsqtot,nlatzones,n1layer,eq_incr)
cTEST
c	print*,ih,"from isqre",i,x(i),y(i)
	   call span(ih,ymi,yma,xmi,xma,nsqrs,nsqtot,nlatzones,eq_incr)
	   ind=(iv-1)*n1layer+ih
cTEST
c	write(*,"(2(i5,2x),4(2x,f12.6))")ih,ind,ymi,yma,xmi,xma
	   if(ind.ne.ind0)then !crossed over to another voxel
	      if(iv.ne.iv0)then !vertical intersection
	            if(z(i).eq.z0)then
	               print*,j,"exception vert"
	               ind=ind0
	               goto11
	            endif
	         ivint=min(iv,iv0)
	         zint=rbnd(ivint)
		 xint=x0+(x(i)-x0)*(zint-z0)/(z(i)-z0)
		 yint=y0+(y(i)-y0)*(zint-z0)/(z(i)-z0)
		 dint=d0+(raydel(i)-d0)*(zint-z0)/(z(i)-z0)
c		 write(*,"(a8,1x,3(f12.6,1x))")"vertical",zint,xint,yint
	      else !horizontal intersection
		     if(dabs(ymi-yma0).lt.toler)then
	            if(y(i).eq.y0)then
	               print*,j,"exception s to n"
	               ind=ind0
	               goto11
	            endif
		    yint=ymi
		    xint=x0+(x(i)-x0)*(yint-y0)/(y(i)-y0)
		    zint=z0+(z(i)-z0)*(yint-y0)/(y(i)-y0)
		    dint=d0+(raydel(i)-d0)*(yint-y0)/(y(i)-y0)
c		    write(*,"(a8,1x,3(f12.6,1x))")"n to s",zint,xint,yint!test
		 elseif(dabs(yma-ymi0).lt.toler)then
	            if(y(i).eq.y0)then
	               print*,j,"exception n to s"
	               ind=ind0
	               goto11
	            endif
		    yint=yma
		    xint=x0+(x(i)-x0)*(yint-y0)/(y(i)-y0)
		    zint=z0+(z(i)-z0)*(yint-y0)/(y(i)-y0)
		    dint=d0+(raydel(i)-d0)*(yint-y0)/(y(i)-y0)
c		    write(*,"(a8,1x,3(f12.6,1x))")"s to n",zint,xint,yint!test
		 elseif(dabs(xmi-xma0).lt.toler)then
	            if(x(i).eq.x0)then
	               print*,j,"exception e to w"
	               ind=ind0
	               goto11
	            endif
	            xint=xmi
		    yint=y0+(y(i)-y0)*(xint-x0)/(x(i)-x0)
		    zint=z0+(z(i)-z0)*(xint-x0)/(x(i)-x0)
		    dint=d0+(raydel(i)-d0)*(xint-x0)/(x(i)-x0)
c		    write(*,"(a8,1x,3(f12.6,1x))")"e to w",zint,xint,yint!test
		 elseif(dabs(xma-xmi0).lt.toler)then
	            if(x(i).eq.x0)then
	               print*,j,"exception w to e"
	               ind=ind0
	               goto11
	            endif
	            xint=xma
		    yint=y0+(y(i)-y0)*(xint-x0)/(x(i)-x0)
		    zint=z0+(z(i)-z0)*(xint-x0)/(x(i)-x0)
		    dint=d0+(raydel(i)-d0)*(xint-x0)/(x(i)-x0)
c		    write(*,"(a8,1x,3(f12.6,1x))")"w to e",zint,xint,yint!test
		 else
		    print*,"problem",j
		    print*,xmi0,xma0,ymi0,yma0
		    print*,xmi,xma,ymi,yma
		    stop "points on ray are too far"
c	            j=j+1
c	            goto10
		 endif
	      endif
cTEST
c	      write(*,"(4(f12.6,2x))")xint,yint,zint,dint
c	      write(99,"(4(f12.6,2x),4(i5,2x))")xint,yint,zint,dint,ih0,iv0,ih,iv

	      zav=(zint+z0)/2.
	      ds=((dint-d0)/radian)*zav
c	      if(ds.ne.0.)then
	      if(ds.ne.0..and.zint.ge.rcmb.and.z0.ge.rcmb)then
	         do l=1,nstartmod-1
	            if(radmod(l).ge.zav.and.zav.ge.radmod(l+1))then
	               a=(vstart(l+1)-vstart(l))/(radmod(l+1)-radmod(l))
	               b= vstart(l)-a*radmod(l)
		       vref=a*zav+b
		       goto33
		    endif
	         enddo
33	         continue

c------------------------------------increment matrix (entries and indices)
	         irec=irec+1
cTEST
c	write(99,"(4(f14.6,1x))")zav,vref,a,b
c	write(89,"(i4,4(f14.6,1x))")l,vstart(l),vstart(l+1),radmod(l),radmod(l+1)
	         write(111,rec=irec)sngl(-ds/vref)
	         write(112,rec=irec)ind0
cTEST
c	write(84,*)ind0,sngl(ds)
	
c	         write(111,*)j,sngl(ds)
c	         write(112,*)ind0
	      endif
	      x0=xint !update x0,y0,z0 
	      y0=yint
	      z0=zint
	      d0=dint
	      iv0=iv !update other variables
	      ih0=ih !not strictly needed
	      ind0=ind 
	      xmi0=xmi
	      xma0=xma
	      ymi0=ymi
	      yma0=yma
	   endif !executed if intersection
11	   continue
	enddo !end of loop over ray path
cTEST
c	pause
c------------------------------------receiver voxel
	if(raydel(npt).ne.d0)then
	   zav=(z(npt)+z0)/2.
	   ds=((raydel(npt)-d0)/radian)*zav
	   irec=irec+1
	   if(ds.ne.0.)then
	      do l=1,nstartmod-1
	         if(radmod(l).ge.zav.and.zav.ge.radmod(l+1))then
	            a=(vstart(l+1)-vstart(l))/(radmod(l+1)-radmod(l))
	            b= vstart(l)-a*radmod(l)
		    vref=a*zav+b
		    goto34
		 endif
	      enddo
34	      continue
cTEST
c	write(99,*)zav,vref
	      write(111,rec=irec)sngl(-ds/vref)
	      write(112,rec=irec)ind
	   endif
cTEST
c	write(84,*)ind0,sngl(ds)
	endif
c------------------------------------increment matrix (pointer) and data vector
	write(113,*)irec
	write(109,*)datum
c--------------------------------------output for synthetic database
        write(50,"(5(f15.8,1x))")elon,elat,depth,slon,slat
	j=j+1
	
c	stop!test

	goto10
20	continue
	end
