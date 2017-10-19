c--does cross-sections through global block models parametrized as 
c--with equal area blocks - like mine.
	parameter(pi=3.1415926536)
	parameter(nodes=180) ! half number of horizontal points for section
	dimension xla(2*nodes),xlo(2*nodes)
	parameter(eq_incr=2.) ! (lateral size of one horizontal pixel)
	parameter(nlatzones=180/eq_incr,refgrid=eq_incr,iswit=0)
c	parameter(nxy=4592) ! number of horizontal pixels (depends on eq_incr)
c	parameter(nxy=6608) ! number of horizontal pixels (depends on eq_incr)
c	parameter(nxy=18332) ! number of horizontal pixels (depends on eq_incr)
c	parameter(nxy=1656) ! number of horizontal pixels (depends on eq_incr)
	parameter(nxy=10316) ! number of horizontal pixels (depends on eq_incr)
	parameter(nlayers=15) ! number of vertical layers
	parameter(ngrid=nxy*nlayers) 
	dimension ave(nlayers)
	real*4 model(ngrid),db(0:nlayers+1),depths(nlayers)
	character name*80,answ*1,answ2*1
	integer nsqrs(nlatzones),nsqtot(nlatzones+1)

	print*,"please specify model file name"
	read(*,"(a80)")name
	ksign=1
	print*,'switch sign of the model? (y/n)'
	read*,answ
	if(answ.eq.'y'.or.answ.eq.'Y')ksign=-1
	print*,'remove average from each layer?'
	read*,answ2

	call param(eq_incr,nsqrs,nsqtot,nlatzones,n1layer,iswit,refgrid)
cTEST
	print*,nsqrs
	print*,n1layer

	open(21,file='cross.dat')  ! output file

c--thickness of layers is constant in this case..
	dpincr=(6371.-3471.)/nlayers
	print*,'constant layer thickness: ',dpincr
	depths(1)=dpincr/2.
	do i=2,nlayers
	   depths(i)=depths(i-1)+dpincr
	enddo

c--define the cross-section
	print*,'pole of the great circle: latitude and longitude(deg)?'
	read*,capth,capph
	capth=90.-capth
	capth=(capth/180.)*pi
	capph=(capph/180.)*pi

c--find a number of equidistant points on a great circle defined by its pole.
c--the formulae employed here can be found by means of a double rotation
c--from the reference frame whose pole is the geographical north to
c--the one whose pole is the pole of the great circle
	open(1,file='gc.dat')
	sf=sin(capph)
	cf=cos(capph)
	ctsf=cos(capth)*sf
	ctcf=cos(capth)*cf
	st=sin(capth)
	rincr=180./nodes
	print*,'vertical increment = ',rincr
	do i=1,nodes
c--by summing 90 to i*rincr in the next line I set the azimuth to 90,
c--i.e. in the output file the 0-point on the great circle corresponds
c--to longitude 90 in the rotated system (the system whose pole is the
c--pole of the great circle)
	   delta=((90.+i*rincr)/180.)*pi
	   theta=abs(acos(-cos(delta)*st))
	   sphi=(ctsf*cos(delta)+cf*sin(delta))/sin(theta)
	   cphi=(ctcf*cos(delta)-sf*sin(delta))/sin(theta)
	   phi=atan2(sphi,cphi)
	   xla(i)=90.-((theta/pi)*180.)
	   xlo(i)=(phi/pi)*180.
	   if(xlo(i).lt.0.)xlo(i)=xlo(i)+360.
	   xla(nodes+i)=-xla(i)
	   xlo(nodes+i)=xlo(i)+180.
	   if(xlo(nodes+i).lt.0.)xlo(nodes+i)=xlo(nodes+i)+360.
	enddo
c--store the great circle in a file
	do i=1,2*nodes
	   write(1,*)xlo(i),xla(i)
	enddo
	close(1)
c------great circle is defined and stored

	do i=1,nlayers
	   ave(i)=0.
	enddo
c--read whole model and store it in memory
	call getmodel(model,name,ngrid,1,nlayers,ksign)
c--compute radii of boundaries between layers
	db(0)=0.
	do idb=1,nlayers
	   db(idb)=db(idb-1)+2.*(depths(idb)-db(idb-1))
	enddo
	do ik=0,nlayers
	   print*,db(ik)
	enddo
c--if required, remove the average from each layer
	if(answ2.eq.'y'.or.answ2.eq.'Y')then
c	   call avera1656(nlayers,ave,name,1)
	   call layav(nlayers,ave,"y",n1layer)
	endif

c-----------------------------------------------------loop over depth 
	incre=50
	ibot=db(nlayers)+0.5
	print*,'inner boundary of mantle is ',ibot,' km deep'
	do ii=0,ibot,incre
c--figure out the layer number
	depth=1.*ii
c	print*,'depth',depth
	if(depth.lt.0.)stop
	if(depth.ge.db(nlayers))then
	   if((depth-db(nlayers)).lt.incre)then
	      depth=db(nlayers)
	      il=nlayers
	      goto221
	   else
	      print*,'too deep!'
	      goto441
	   endif
	endif
	do idb=1,nlayers
	if(depth.ge.db(idb-1).and.depth.lt.db(idb))then
	   il=idb
	   goto221
	endif
	enddo
 221	continue
c	print*,'layer n. ',il,', increment=',rincr
c--for this depth, evaluate the model along the whole great circle
	do i=1,nodes*2
	   xlon=xlo(i)
	   colat=90.-xla(i)
	   call readmodel(model,xlon,colat,il,value,
     &	                  n1layer,ngrid,nsqrs,nsqtot,nlatzones,eq_incr)
	   write(21,*)i*rincr,6371.-depth,value-(ave(il)*ksign)
	enddo
	enddo
c----------------------------------------------end of loop over depth

441	close(21)
	end


C..................................................................
C..................................................................
	subroutine getmodel(model,name,n,io,nlayers,ksign)
c--reads a block model from file and stores it in memory.
c--the parametrization of the model is describes by input variables.
	character name*80
	real*4 model(n)

	print*,'sign is',ksign
	print*,"model has",n," coefficients"
	print*,'reading from model',name

	open(unit=io,file=name,status='old')

	read(io,*)
	do i=1,n
	   read(io,*)k,model(i)
	   model(i)=ksign*model(i)
	   if(i.ne.k)stop "mismatch in coefficient index"
	enddo
	close(io)
	return
	end
c==================================================================
	subroutine readmodel(model,lon,colat,layern,value,
     &	                     nxy,n,nsqrs,nsqtot,nlatzones,eq_incr)
	real*4 model(n),lon
	real*8 dlat,dlon
	dimension nsqrs(nlatzones),nsqtot(nlatzones+1)

	ntot=nxy*(layern-1)
	dlat=90-colat
	dlon=lon
	ntot=ntot+isqre(dlat,dlon,nsqrs,nsqtot,nlatzones,n,eq_incr)
	value=model(ntot)
	return
	end
c==================================================================
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
	return
	end
c==================================================================
	subroutine layav(nlayers,ave,answl,n1layer)
	dimension ave(nlayers)
	character*1 answl
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
	enddo
	print*,'I computed all the ',nlayers,' averages'
	rewind(1)
	return
	end
