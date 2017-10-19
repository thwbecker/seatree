	parameter(pi=3.141592653589)
	character cbin_mid*3,chataid*7,chataid2*7,rayin*2
c--by setting iswiresc to 0 we apply no rescaling...
	parameter(iswiresc=0)
c--set isotro=1 : isotropic; isotro=2 : radially anisotropic
	parameter(isotro=1)
c--if delay time > cutoff then assign a smaller weight.
	parameter(cutoff=3.)
	real*4 datum
	parameter(nlay=15)
	character*256 raynam,filnam
	real*8 ddel,ddep,dpres1
	real*8 rboundary(0:20)
	real*4 rescale(20)

	parameter(iswit=0)
	parameter(nlatzomax=180)
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
	common/deltastats/deltami,deltama

	deltami=180.
	deltama=0.

	print*,"size of pixels?"
	read*,eq_incr
	refgrid=eq_incr*1.
	nlatzones=180./eq_incr
	if(nlatzones.gt.nlatzomax)stop "coarse pixels too small"

	print*,"start from datum n.?"
	read*,jfirst
	print*,"end at datum n.?"
	read*,jlast
        write(chataid,'(i7.7)')jfirst
        write(chataid2,'(i7.7)')jlast
        print*,"what data? (P0,bc,df,ab,cp)"
        read*,rayin
	print*,'what bin?'
	print*,'min. delta (degrees)?'
	read*,bin_min
	print*,'max. delta (degrees)?'
	read*,bin_max
	write(cbin_mid,'(i3.3)')int((bin_min+bin_max)/2.)
	write(*,*)"bin:",bin_min,bin_max
	print*,'cutoff=',cutoff

	geoco=0.993277
	radian=180./pi
	ncoltotal=0
	nto1=0
	nto2=0
	nto3=0
	nto4=0
	nto5=0
	nto6=0
	nto7=0
	nto8=0
	nto9=0
	nto10=0
	ndata=0
	iflapkp=0
	kflag=0
	iftstcor=0

c-----define parameterization
c	call param(eq_incr,nsqrs,nsqtot,nlatzones,n1layer,iswit,refgrid)
	call param(eq_incr,nsqrs,nsqtot,nlatzones,n1layer,iswit,refgrid
     &,nlatzomax)

c	open(94,file='check.d')

	open(74,file='core.test')
c--name of the reference model (anisotropic prem) file
	filnam='aniprmc'

	print*,"opening data files"
	print*,'../hrvdata/data.'//rayin//'.smm'

c--select ray
        open(33,file='../hrvdata/data.'//rayin//'.smm.ascii',status="old")
c     &  access='direct',recl=4,status='old',form='unformatted')
        open(34,file='../hrvdata/quakes.'//rayin//'.smm.ascii',status="old")
c     &  access='direct',recl=12,status='old',form='unformatted')
        open(35,file='../hrvdata/receiv.'//rayin//'.smm.ascii',status="old")
c     &  access='direct',recl=8,status='old',form='unformatted')

        iflapkp=0
        if(rayin.eq."P0")then
           raynam="P"
        elseif(rayin.eq."bc")then
           raynam="PKP"
           iflapkp=1
        elseif(rayin.eq."ab")then
           raynam="PKP"
           iflapkp=2
        elseif(rayin.eq."cp")then
           raynam="PcP"
        elseif(rayin.eq."df")then
           raynam="PKIKP"
        else
           stop"cant deal with this phase"
        endif

cTEST
	print*,"opening matrix files"
c--open the files where the MATRIX will be stored
	raynam='P'
	open(11,file='a.vx_mat_'//rayin//'_'//chataid//'.'//chataid2,
     &  access='direct',recl=4,form='unformatted')
        open(12,file='a.vx_ind_'//rayin//'_'//chataid//'.'//chataid2,
     &  access='direct',recl=4,form='unformatted')
        open(13,file='a.vx_pnt_'//rayin//'_'//chataid//'.'//chataid2)
        open(9,file='d.vx_vec_'//rayin//'_'//chataid//'.'//chataid2)

c--read summary data
	j=0
	do j=jfirst,jlast
cTEST
c	print*,"record ",j
c-----------this version will crash if end-of-file is hit
	read(33,*)datum
	read(34,*)sum_soulon,sum_soulat,sum_soudep
	read(35,*)sum_stalon,sum_stalat
cTEST
c	print*,sum_soulon,sum_soulat,sum_soudep,sum_stalon,sum_stalat

	if(sum_soulon.lt.0.)sum_soulon=sum_soulon+360.
	if(sum_stalon.lt.0.)sum_stalon=sum_stalon+360.
	elon=sum_soulon
	elat=sum_soulat
	depth=sum_soudep
	slon=sum_stalon
	slat=sum_stalat


cTEST
c	print*,"saving source station locations"
c	open(884,file="sost.txt")
c	write(884,*)slon,slat
c	write(884,*)elon,elat
c	close(884)

c--call all the programs that build up the matrix

c--here latitudes are defined correctly as latitudes and NOT co-latitudes.
c--longitude is defined so that it takes also negative values.
c	call delaz1(elat,elon,slat,slon,delta,azep,azst)
c--delazs IS A VERSION OF SAME ROUTINE I'VE GOTTEN FROM JOHN W. IN 2003.
c--I am not sure it is the same as delaz1 (what about geocentric correction),
c--this will need some tests!
	call delazs(elat,elon,slat,slon,delta,azep,azst)
cTEST
c	print*,elat,elon,slat,slon
c	print*,"epicentral distance=",delta
c	pause

	ndata=ndata+1
c--check if datum is to be discarded
	if(abs(pres1).gt.10.)then
	   nto1=nto1+1
	   goto10
	endif
	if(delta.le.25.)then
cTEST
c	print*,"small delta=",delta
	   nto6=nto6+1
	   goto10
	endif

	if((delta.gt.bin_max).OR.(delta.le.bin_min))then
	   nto10=nto10+1
	   goto10
	endif

	ddel=delta*1.d0
	ddep=depth*1.d0
	nto6old=nto6
	nto7old=nto7

	call find_path(ddel,raynam,ddep,1,0,filnam,nto6,nto7,iflapkp)
	if(nto6.gt.nto6old)goto10
	if(nto7.gt.nto7old)goto10
c--the input parameters mean respectively: delta; name of phase;
c--depth of source; a number that has to be non zero if want derivatives;
c--a number that is =0 for anisotropic layer, =1 for equiv. isotropic model
c--(PREM); name of such reference model (anisotropic PREM in general).
c--nto6 and nto7 are the number of data thrown out because of
c--different malfunctions of ttpath (output).

c--define radial parametrization (only once)
	if(kflag.eq.0)then
	   call parametriz(nlay,rboundary,rescale,iswiresc)
	   kflag=1
	endif

c--for rays that travel through the inner core, correct
c--the datum for the effect of inner core anisotropy.
	if(raynam.eq.'PKIKP')then
c--I will assume that Weijia's subroutine wants geocentric coordinates.
c--Indeed since we are using summary rays it should not make a big
c--difference...however:
           slat1=radian*ATAN(GEOCO*TAN(dble(sLAT)/radian))
	   slon1=slon
	   if(slon.gt.180.)slon1=slon-360.
           elat1=radian*ATAN(GEOCO*TAN(dble(ELAT)/radian))
	   elon1=elon
	   if(elon.gt.180.)elon1=elon-360.
	   CALL aniso_cor0(SLAT1,SLON1,0.,ELAT1,ELON1,depth,RESIC,IERR)
	   datum=datum-RESIC

c-------------------------------------------------------
c--check the inner core anisotropy as given by aniso_cor
	   IF(iftstcor.EQ.0)THEN
	      iftstcor=1
	      TSTLAT1=SLAT1
	      TSTLON1=SLON1
	      WRITE(74,*)TSTLON1,TSTLAT1
	   ENDIF

	diffla=abs(slat1-tstlat1)
	difflo=abs(slon1-tstlon1)
	   if((diffla.le.5.).and.(difflo.le.5.))then
	      WRITE(74,*)eLON1,eLAT1,RESIC
	   endif
c-------------------------------------------------------

	endif

	dpres1=datum*1.d0
	deltamiold=deltami
	deltamaold=deltama
	call path_3d(n1layer,nsqrs,nsqtot,nlatzones,
     &ddep,ddel,ELAT,ELON,SLAT,SLON,dpres1,rboundary,
     &rescale,nlay,ncoltotal,nto2,nto3,nto4,nto5,nto8,nto9,
     &cutoff,isotro,eq_incr,nlatzomax)
c	if(deltamiold.ne.deltami)print*,"min delta=",deltami
c	if(deltamaold.ne.deltama)print*,"max delta=",deltama

10	continue
	if(mod(ndata,1000).eq.0)then
c	   open(44,file='stats.d'//cbin_mid)
	   open(44,file='stats.d')
	   write(44,*)"bin:",bin_min,bin_max
	   write(44,*)"raw number of data:",ndata
	   write(44,*)'number of data thrown out...'
	   write(44,*)'...because of too large delay time',nto1
	   write(44,*)'...problems with interpolation',nto2
	   write(44,*)'...screwed up path',nto3
	   write(44,*)'...repeated indexes',nto4
	   write(44,*)'...delta-finding procedure did not converge',nto5
	   write(44,*)'...delta outside range',nto6+nto10,nto6,nto10
	   write(44,*)'...zero-finding did not converge (ttpath)',nto7
	   write(44,*)'...problems in kernels integration',nto8
	   write(44,*)'...delta too large; equations break down',nto9
	   ntrounaut=nto1+nto2+nto3+nto4+nto5+nto6+nto7+nto8+nto9+nto10
	   write(44,*)'number of data used:',ndata-ntrounaut
	   close(44)
	endif

c-----------------------------------------------------END OF THE MAIN LOOP
	enddo
105	print*,'total number of couples:',j-jfirst
c--final statistics:
c	open(45,file='stats.d'//cbin_mid)
	open(45,file='stats.d')
	write(45,*)"bin:",bin_min,bin_max
	write(45,*)"raw number of data:",ndata
	write(45,*)'number of data thrown out...'
	write(45,*)'...because of too large delay time',nto1
	write(45,*)'...problems with interpolation',nto2
	write(45,*)'...screwed up path',nto3
	write(45,*)'...repeated indexes',nto4
	write(45,*)'...delta-finding procedure did not converge',nto5
	write(45,*)'...delta outside range',nto6+nto10,nto6,nto10
	write(45,*)'...zero-finding did not converge (ttpath)',nto7
	write(45,*)'...delta too large; equations break down',nto9
	ntrounaut=nto1+nto2+nto3+nto4+nto5+nto6+nto7+nto8+nto9+nto10
	write(45,*)'number of data used:',ndata-ntrounaut
	close(45)
	close(1)
	close(2)
	close(33)
c	close(922)
c	close(923)
	close(9)
	close(11)
	close(411)
	close(12)
	close(13)
	open(25,file="num.obs")
	write(25,*)ndata-ntrounaut
	close(25)
	end

	subroutine coordsuper(NBLOC,blocla,bloclo,nsqrs,nlatzones,eq_incr)
c--given a cell index on the Earth's surface, finds longitude and latitude
c--(NOT colatitude) of its center.
	dimension nsqrs(nlatzones)
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
