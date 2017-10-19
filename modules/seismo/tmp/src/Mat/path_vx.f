c    this program computes travel times,etc and
c derivatives of travel time w.r.t. the model parameters (for a
c prem-like model parametrization currently).
c
c  this is a major rewrite of tpg50 by guy masters, oct 97
c  i tried to remove the model from the body of the code but this
c  is not quite perfect. I added readmod to read the model, getmod
c  to evaluate the model at a specific radius, and cder which
c  evaluates the derivative of a model parameter. the place
c  where the model parameterization remains is in fqanis which
c  makes assumptions about the parameterization when integrating
c  the model derivatives
c
c  i also added raynam,ichdec,addray which allows the user to
c  input standard ray names -- the user is prompted for SH
c  or SV if the ray has only S legs
c
c  all of the remaining routines have been rewritten to some degree
c  even though their names are the same as before. Some have been
c  modified extensively but most changes are to make the code more
c  linear
c
c  Results are written to file ttnew.out and some are printed to
c  the screen. The path is written to file "path" and is in x,y
c  format where y=rcos(delta) and x=rsin(delta) and r is normalized
c  so that the surface of the earth is r=1. "path" also includes
C  3 circles of radius ricb, rmcb, and 1 respectively -- each
c  circle comprises 361 points so the first 1083 points of the
c  file are the circles. Each ray path is separated by x,y=9999,9999
c  (Subroutines raypath and addpath have been added to do the path
c  computation)
c  
c
c  I have changed the model format to be consistent with the
c  mode code
c
c  I added the computation of t*
c
c    the user inputs:
c dep    - source depth in km
c ray name (will ask of SH or SV for pure S rays)
c ider   - =1 if derivatives are to be calculated
c iso    - =0(isotropic layers remain isotropic)/ =1(otherwise)
c imod   - =0(anisotropic model is used)/ =1(averaged isotropic model)
c iopt   - =0 for a table with evenly spaced deltas
c          =1 for a particular ray parameter p
c          =2 for a particular delta
c          =3 for a table with evenly spaced p's
c  iopt=0 or 3 will cause a prompt for table spacing
c

	SUBROUTINE FIND_PATH(DELTA,RAY,DEP,IDERSV,Jiso,filnam,nto6,nto7,
     &iflapkp,ipolar)

      implicit double precision (a-h,o-z)
      character*256 ray,filnam
      dimension qvec(500)
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/slow/picblp,picbls,picbup,pcmblp,pcmbup,pcmbus,
     +   p660lp,p660ls,lic,loc,llm,lsrc,klay(20)
      common/ray$/pmin,pmax,nps(20,2),npsrc(2),nn,iup,jref,npleg,nsleg
      external f

c------------------------------Lapo(11/98)------------------------------
	COMMON/CMBKER/cmbkernl
	data pi/3.14159265358979d0/
	common/SLOWLAPO/pcmblpLAPO,pcmbupLAPO
c------------------------------Lapo(11/98)------------------------------

	common/checkzero/ntoerr

c      open(9,file='ttnew.out')
c      open(10,file='path')

c---------------------------------LAPO----------------------------------
c      open(121,file='time.sub')
c      open(125,file='radi.sub')
c	open(124,file='tzech')
c---------------------------------LAPO----------------------------------

c
c***** set the following without asking
c  order of GL integration
      idleg=6
c  iso=0 means we compute isotropic derivs for isotropic prem-layers 
      iso=1
c imod =0 for anisotropic layer, =1 for equiv. isotropic model
C      imod=0
	IMOD=JISO
c************************************
c
      call readmod(filnam)

c--now pcmblp,pcmbup are the reference P-slownesses just below and just 
c--above the cmb, respectively.
c--xt(loc) is the cmb radius.

c      print *,'input source depth in km'
c      read(*,*) dep
c*** we must set up source layer before calling raynam
      lsrc=0
      rdep=6371.d0-dep
      do i=1,numlyr
        if(xb(i).lt.rdep.and.xt(i).ge.rdep) lsrc=i
      enddo
C      ray=' '
C      print *,'enter ray name'
C      read(*,*) ray
      call raynam(ray,ierr)
      if(ierr.ne.0) then
        print *,'error in ray name'
        stop
      end if
      ish=0

      if(ifanis.ne.0) then
c*** check to see if this is pure S -- then ask if SH or SV
        if(npleg.eq.0) then
c           print *,'enter 1 for SH or 0 for SV'
c           read(*,*) ish
	ish=ipolar
           if(ish.ne.1) ish=0
        end if
      end if

c*** ask for these control variables 
C      print *,'enter ider (non zero if want derivatives): '
C      read(*,*) idersv
C      print *,'enter iopt (=0, table of delta, =1 a single p'
C      print *,'=2 a single delta, =3 a table of ps): '
C      read(*,*) iopt
C--SET IOPT WITHOUT ASKING (A SINGLE DELTA)
	IOPT=2


      if(iopt.eq.1) goto 140
      if(iopt.eq.2) goto 150
      if(iopt.eq.3) goto 160
c ===========================================================
c ----iopt=0, i.e. we asked for a table with evenly spaced deltas
  120 p1=pmin
      p2=pmax
c ----in d1 and d2 we put the values of d(delta)/d(p) for p=p1 and
c ----p=p2; if they have opposite signs, d(delta)/d(p) must be zero in
c ----between. the value for which this happens is found by a call
c ----to subroutine zero then the user is prompted for a new p range
      itim=2
      ider=0
      call deriv(p1,qvec,xtu)
      delex1=qvec(1)
      d1=qvec(2)
      call deriv(p2,qvec,xtu)
      delex2=qvec(1)
      d2=qvec(2)
      if(d1*d2.lt.0.d0) then
        ptarg=0.d0
        call zero(p3,p1,p2,f,ptarg)
        print *,'old pmin and pmax : ',pmin,pmax
        print *,'d(delta)/d(p) vanishes at p= ',p3
        print *,'enter new pmin,pmax (neg vals to quit)'
        read(*,*) pmin,pmax
        if(pmin.lt.0.d0.or.pmax.lt.0.d0) stop
        goto 120
      end if
c*** monotonic
      print *,'input spacing of table by delta in degs'
      read(*,*) delf
      if(delf.le.0.d0) goto 1000
      dismin=delex1
      dismax=delex2
      if(delex1.gt.delex2) then
        p1=pmax
        p2=pmin
        dismin=delex2
        dismax=delex1
      end if
      delreq=-delf
  130   delreq=delreq+delf
        if(delreq.lt.dismin) goto 130 
        if(delreq.gt.dismax) goto 1000
        itim=1
        ider=0
        call zero(pp,p1,p2,f,delreq)
        if(pp.lt.0.d0) goto 130
        itim=3
        ider=idersv
        call deriv(pp,qvec,xtu)
        call prray(ray,pp,qvec(1),qvec(5),xtu)
        call raypath(ray)
        p1=pp
      goto 130
c ===================================================================
c this block is called if we want particular values of p
  140   print *,'suggested p range for this phase :  ',pmin,pmax
  145   print *,'enter p (neg to quit): '
        read(*,*) pp
        if(pp.lt.0.d0) goto 1000
        itim=3
        ider=idersv
        call deriv(pp,qvec,xtu)
        call prray(ray,pp,qvec(1),qvec(5),xtu)
        call raypath(ray)
      goto 145
c ==================================================================
c--THE FOLLOWING BLOCK IS WHAT I WANT!
c this block converts a request for a value of delta into a request for
c a value of a specific p
  150   p1=pmin
        p2=pmax
        itim=2
        ider=0
        call deriv(p1,qvec,xtu)
        delex1=qvec(1)
        d1=qvec(2)
        call deriv(p2,qvec,xtu)
        delex2=qvec(1)
        d2=qvec(2)
        if(d1*d2.lt.0.d0) then
          ptarg=0.d0
          call zero(p3,p1,p2,f,ptarg)
c          print *,'old pmin and pmax : ',pmin,pmax
c          print *,'d(delta)/d(p) vanishes at p= ',p3
c          print *,'enter new pmin,pmax (neg vals to quit)'
c          read(*,*) pmin,pmax

c---------------------------------LAPO----------------------------------
c	   print*,pmin,pmax
c--automatically reset pmax and pmin for PKP bc and df.
c--PKPbc:
	   if(iflapkp.eq.1)pmax=p3-1.d-4

c--PKPab
	   if(iflapkp.eq.2)pmin=p3+1.d-4
c	   print*,pmin,pmax
c---------------------------------LAPO----------------------------------

          if(pmin.lt.0.d0.or.pmax.lt.0.d0) stop
          goto 150
        end if

c---------------------------------old-----------------------------------
c        print *,"delta ranges from ",delex2," to ",delex1
c        print *,"as p decreases from ",p2," to ",p1
c  155   print *,"type in delta (negative to quit)"
c        read *,delreq
c---------------------------------old-----------------------------------
c---------------------------------LAPO----------------------------------
	delreq=delta
c	print*,delex2,delta,delex1,DEP

	if((iflapkp.eq.2).and.
     &((delta.lt.delex1).or.(delta.gt.delex2)))then
	   print*,'no!'
	   nto6=nto6+1
	   goto1000
	elseif((iflapkp.eq.1).and.
     &((delta.lt.delex2).or.(delta.gt.delex1)))then
	   print*,'no!'
	   nto6=nto6+1
	   goto1000
	endif
c---------------------------------LAPO----------------------------------

        if(delreq.lt.0.d0) go to 1000
        itim=1
        ider=0
        call zero(pp,p1,p2,f,delreq)

c---------------------------------LAPO----------------------------------
	if(ntoerr.eq.1)then
	   nto7=nto7+1
	   return
	endif
c---------------------------------LAPO----------------------------------

        itim=3
        ider=idersv
        call deriv(pp,qvec,xtu)
        call prray(ray,pp,qvec(1),qvec(5),xtu)

c------------------------------Lapo(11/98)------------------------------
c--at this point we know ray parameter, turning point, reference slowness 
c--just below and just above cmb, and reference cmb radius.
c--pp is in seconds/degree
	pplapo=(180.d0/pi)*pp
c	print*,'turning point radius:',xtu
c	print*,'ray parameter:',pp,xt(loc)/pp
c	print*,'ray parameter, rad:',pplapo,xt(loc)/pplapo
c	print*,'slowness above and below cmb:',pcmbuplapo,pcmblplapo
c	print*,'CHECK THIS',
c     &(180.d0/pi)*Asin(pplapo*(1.d0/pcmbuplapo)/xt(loc))
c	print*,'vel. above and below cmb:',1.d0/pcmbuplapo,
c     &1.d0/pcmblplapo
c	print*,'cmb radius:',xt(loc)
	if(xtu.lt.xt(loc))then
c--assign a value to cmb structure kernel (Morelli & Dziewonski, 1987).
c--I want the elements of the A matrix to have all the same
c--dimension (seconds), so I define the kernel consistently, and
c--my unknown will be the RELATIVE cmb-perturbation (deltar/r)
c	   cmbtauabo=pplapo*pplapo*(1.d0/(xt(loc)*xt(loc)))
c	   cmbtauabo=(pcmbuplapo*pcmbuplapo)-cmbtauabo
	   cmbtauabo=pplapo*pplapo
          cmbtauabo=(pcmbuplapo*pcmbuplapo*xt(loc)*xt(loc))
     &-cmbtauabo

	   cmbtauabo=sqrt(cmbtauabo)

c	   cmbtaubel=pplapo*pplapo*(1.d0/(xt(loc)*xt(loc)))
c	   cmbtaubel=(pcmblplapo*pcmblplapo)-cmbtaubel
	   cmbtaubel=pplapo*pplapo
	   cmbtaubel=(pcmblplapo*pcmblplapo*xt(loc)*xt(loc))
     &-cmbtaubel

	   cmbtaubel=sqrt(cmbtaubel)

	   cmbkernl=-(cmbtauabo-cmbtaubel)

c	   print*,'cmb structure kernel',cmbkernl
c	   pause
	else
	   cmbkernl=0.d0
	endif
c	pause
c------------------------------Lapo(11/98)------------------------------

        call raypath(ray)
c	close(10)
c---------------------------------old-----------------------------------
c      goto 155
c---------------------------------old-----------------------------------
c---------------------------------LAPO----------------------------------
	goto1000
c---------------------------------LAPO----------------------------------

c =======================================================================
c this block is called if we want a table of p's
  160   print *,'suggested p range for this phase :  ',pmin,pmax
        print *,'enter pmin,pmax, and number of p values : '
        read(*,*) p1,p2,npx
        if(npx.gt.1) then
          dpx=(p2-p1)/(npx-1)
          do i=1,npx
            pp=p1+(i-1)*dpx
            itim=3
            ider=idersv
            call deriv(pp,qvec,xtu)
            call prray(ray,pp,qvec(1),qvec(5),xtu)
            call raypath(ray)
          enddo
        end if
c ==================================================================
 1000 continue

c	close(125)
c	close(9)
c	close(124)
      RETURN

      end

c=================================================================
      subroutine prray(ray,pp,qvec,qray,xtu)
c*** prints out results for this ray type and writes to unit 9
c*** qvec(1)=delta, qvec(2)=dX/dp, qvec(3)=time, qvec(4)=tstar
      implicit real*8(a-h,o-z)
      character*(*) ray
      dimension qray(4,5,20),qvec(*)

c---------------------------------LAPO----------------------------------
	COMMON/TEMPO/PREMTIME
c---------------------------------LAPO----------------------------------

      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/slow/picblp,picbls,picbup,pcmblp,pcmbup,pcmbus,
     +   p660lp,p660ls,lic,loc,llm,lsrc,klay(20)
      data if1/1/
      if(if1.eq.1) then
c        print 930
c        write(9,930)
  930   format(3x,'Ray',8x,'Distance',5x,'Time',10x,'t*',12x,'dXdp',9x,
     +   'p',5x,'Bot. rad')
        if1=0
      end if
      qdel=qvec(1)
      qddp=qvec(2)
      qtim=qvec(3)
      qtst=qvec(4)

c---------------------------------LAPO----------------------------------
	premtime=qtim
c---------------------------------LAPO----------------------------------

c        write(9,930)
c      write(9,920)ray(1:12),qdel,qtim,qtst,qddp,pp,xtu
c      print 920,ray(1:12),qdel,qtim,qtst,qddp,pp,xtu
      if(ider.ne.0)then
        do 170 l=1,numlyr
c          write(9,900)l 
          do 170 k=1,4
c            write(9,910) k,(qray(k,j,l),j=1,5)
170	  continue
      end if
  900 format(3x,28("-"),"   l a y e r   ",i2,3x,28("-")
     +  /,3x,"k",8x,"vph",12x,"vpv",12x,"vsh",12x,"vsv",12x,"eta")
  910 format(2x,i2,5(3x,e13.6))
  920 format(1x,a12,3(f8.3,4x),f12.4,4x,f7.4,f10.3)
      return
      end

      subroutine deriv(pp,qvec,xtu)
c     deriv accepts as input the
c     ray parameter pp in seconds/degree and returns array qvec:
c   qvec(1)=distance; qvec(2)=dxdp; qvec(3)=time; qvec(4)= tstar
c   ...rest are derivatives of travel time w.r.t. the model 
c   coefficients for each of vph,vpv,vsh,vsv,and eta - there are
C   5 physical params each described by a cubic (for prem) so there 
C   are 20(=nplay) model parameters per layer
c   xtu is the turning point radius 
c
c     further input must be passed by the calling routine through
c     common blocks:
c         common/layr/nl(20),xb(20),xt(20),ifanis,nplay
c         common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
c     where
c         nl(i) = number of levels in layer i, xb(i) and xt(i) are
c                 bottom and top in kms. of layer i
c         rdep=source radius in kilometers
c         itype=1,2,3 for sh,p,sv waves
c         imod=0 for anisotropic model,1 for isotropic
c         itim: =1 to calculate only the deltas
c               =2 to calculate deltas and dXdp
c               =3 to calculate deltas, travel times, and deriva-
c                  tives of delta w.r.t. p and tstar
c          ider=1 to calculate model derivs
c         idleg=degree of the legendre polynomials to be used in
c                 the integrations
c         rnorm=length scale with respect to which the polynomial
c                 coefficients are normalized
c         numlyr=number of layers in the model - up to 20
c
      implicit double precision (a-h,o-z)
      logical isplit
c*** qints must be dimensioned large enough to be nplay+4
      dimension qvec(*),qints(24)
      common/path$/delt(2,800,2),nfin(2),kdep
      common/ray$/pmin,pmax,nps(20,2),npsrc(2),nn,iup,jref,npleg,nsleg
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/slow/picblp,picbls,picbup,pcmblp,pcmbup,pcmbus,
     +   p660lp,p660ls,lic,loc,llm,lsrc,klay(20)
      common/xnorm/xnorm
      data pi/3.14159265358979d0/
      rad=180.d0/pi
      p=pp*rad
      npar=numlyr*nplay+4
      do j=1,npar
        qvec(j)=0.d0
      enddo
c*** this is loop over ray type
      do 1000 ijk=1,2
      if(ijk.eq.1.and.npleg.eq.0) goto 1000
      if(ijk.eq.2.and.nsleg.eq.0) goto 1000
c*** itype=1 for SH, =2 for P, = 3 for SV
      itype=ijk+1
      if(ijk.eq.2.and.ish.ne.0) itype=1
c*** jpath counts increments in delta used in making path
      jpath=0
c ===================================================================
c ----loop over all the layers crossed by the ray
c isplit(logical): true if the source is in the current (n-th) shell
c rdep  : distance of the epicenter from the center of the earth
c iflag : =0(if ray traversed entire layer n); =1(if ray is reflected
c         at the top of the layer n); =2(if ray turned within layer n);
c         its value returns from the call to subroutine integr
      do 100 n=1,nn
        nlast=n
        if(n.ne.1) jpath=klay(n-1)
        if(nps(n,ijk).eq.0) go to 100
        isplit=.false.
        if(n.eq.lsrc) then
          xbsav=xb(n)
          xtsav=xt(n)
          if(rdep.ne.xt(n))then
c    if the focus is in this shell, the bottom of the shell is
c    set equal to rdep and the upper half-shell is done, the
c    lower half-shell is done by means of the goto 105 at the
C    end of the do 100 block
             xb(n)=rdep
             isplit=.true.

c---------------------------------LAPO----------------------------------
c	   write(124,*)'layer just above the source:'
c---------------------------------LAPO----------------------------------

          endif
        end if
  105   continue

c---------------------------------LAPO----------------------------------
c	write(124,*)'call integr',klay(n-1),jpath
c---------------------------------LAPO----------------------------------

        call integr(n,qints,iflag,jpath,ijk)
c -------ray reflected at the top of layer n
        if(iflag.eq.1) then
           nlast=n-1
           xtu=xt(n)
           go to 1000
        endif
        fmult=nps(n,ijk)
        if(n.eq.lsrc.and.isplit) fmult=npsrc(ijk)
        do i=1,4
          qvec(i)=qvec(i)+qints(i)*fmult
        enddo

c---------------------------------LAPO----------------------------------
c	write(124,*)'fmult=',fmult
c---------------------------------LAPO---------------------------------- 

        if(ider.ne.0) then
          ind=(n-1)*nplay+4
          do i=1,nplay
            qvec(i+ind)=qvec(i+ind)+qints(4+i)*fmult
          enddo
        end if
c*** if this is split source layer -- go back for bottom bit
        if(n.eq.lsrc) then
          xb(n)=xbsav
          if(isplit)then
            xt(n)=rdep
            isplit=.false.

c---------------------------------LAPO----------------------------------
c	   write(124,*)'layer just below the source:'
c---------------------------------LAPO---------------------------------- 

            goto 105 
          else
            xt(n)=xtsav
          endif
       end if
c ---- if flag=2 ray turned in this layer; xtu (turning point
c ----radius) is set equal to xnorm (through common block /xnorm/)
        if(iflag.eq.2)then
          xtu=xnorm
          goto 1000
        endif
  100   continue
c ======================================================================
 1000 nfin(ijk)=jpath
c*** end of loop over wave type
      qvec(1)=qvec(1)*rad
      qvec(2)=qvec(2)*rad*rad
      if(jref.eq.1) xtu=xt(loc)
      if(jref.eq.2) xtu=xt(lic)
      return
      end

      subroutine integr(n,qints,iflag,jpath,i12)
c   subroutine to integrate contributions to delta,time, and
c   time derivatives (returned in qints) for layer n.
c     qints(1)=d
c     qints(2)=ddp
c     qints(3)=t
c     qints(4)=tstar
c     qints(5...)=tt model derivatives
c     iflag  =0(if ray traversed entire layer)
c            =1(if ray was reflected at the top of the layer)
c            =2(if ray turned in this layer)
c
      implicit double precision (a-h,o-z)
      dimension qints(*),ind(3)

c---------------------------------LAPO----------------------------------
	COMMON/TIMEINCR/TIMEARRAY(800,2)
c------------------------------LAPO(10/98)------------------------------
c	dimension aktemp(5)
	COMMON/ANKERINC/DQDVPH(800,2),DQDVPV(800,2),DQDVSH(800,2),
     &DQDVSV(800,2),DQDETA(800,2)
C	common/lapoanis/onewayker1(800),onewayker2(800),onewayker3(800),
C     &onewayker4(800),onewayker5(800),rayker(800)
c------------------------------LAPO(10/98)------------------------------
c---------------------------------LAPO----------------------------------

      common/path$/delt(2,800,2),nfin(2),kdep
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/parm1/xa,xc,xn,xl,xf,tdif(5)
      common/xnorm/xnorm
      common/dx/dxtdp,fxtpp,ifl
      external fqanis
      data ind/1,2,4/
      data pi/3.14159265358979d0/
c--num =1(only deltas); =2(delta and dX/dp); =4 (+time and tstar)
      num=ind(itim)
      if(ider.ne.0) num=4+nplay
      do i=1,num
        qints(i)=0.d0
      enddo
      if(rdep.eq.6371.d0) kdep=1
      ifl=-1
      iflag=0
      call qtau(xt(n),n,q1)
      if(q1.lt.0.d0)then 
c*** oops, ray already turned
        iflag=1
        return
      endif
      call qtau(xb(n),n,q1)
c*** find turning point -- assumes 1 per layer!
      if(q1.le.0.d0) then
        call zero1(xnorm,xt(n),xb(n),n,ierr)
        ifl=0
        iflag=2
c ----dxt returns the derivative (dxtdp) of the turning radius w.r.t.
c ----parameter p. for fxtpp, see woodhouse & girnius, iii.14 
        if(itim.ge.2) call dxt(xnorm,n,dxtdp,fxtpp)
      else
        xnorm=xb(n)
      end if
      nlev=nl(n)-1
      xdec=(xt(n)-xb(n))/nlev
      xdecs=xdec+1.
      x1=xt(n)-xnorm
      sumd=0.

c---------------------------------LAPO----------------------------------
	sumti=0.
c------------------------------LAPO(10/98)------------------------------
	sumdqdvph=0.
	sumdqdvpv=0.
	sumdqdvsh=0.
	sumdqdvsv=0.
	sumdqdeta=0.
c------------------------------LAPO(10/98)------------------------------
c---------------------------------LAPO----------------------------------

  10    x2=x1
c*** need fine spacing near center for deep turning rays
        if(xnorm.le.10.d0.and.x2.lt.xdecs) xdec=xnorm
        x1=x1-xdec
        if(x2.le.0.d0) go to 20
        if(x1.le.0.d0) x1=0.d0
ccccccccccccccccccccccccccccccccccccccccc
c**** temporary code to write kernels to unit 99
        if(ider.ne.0) then
          y=0.5d0*(x1+x2)+xnorm
          call qtau(y,n,q1)
          qq=sqrt(q1)
          call tder(y,qq)
c          write(99,900) y,(tdif(k),k=1,5)
  900     format(6g14.6)

c---------------------------------LAPO----------------------------------
c------------------------------LAPO(10/98)------------------------------
c	  do ktemp=1,5
c	     aktemp(ktemp)=tdif(ktemp)
c	  enddo
c------------------------------LAPO(10/98)------------------------------
c---------------------------------LAPO----------------------------------

        end if
cccccccccccccccccccccccccccccccccccccccccc
        xx1=dsqrt(x1)
        xx2=dsqrt(x2)
c ----gauslv is the routine which performs the gauss-legendre integration.
        call gauslv(xx1,xx2,idleg,fqanis,n,qints,num)
        if(xdec.eq.xnorm.and.x1.ne.0.d0) goto 10
c*** save increments in delta as a function of radius
        jpath=jpath+1
        delt(1,jpath,i12)=x1+xnorm
        delt(2,jpath,i12)=qints(1)-sumd

c---------------------------------LAPO----------------------------------
c--save anisotropic kernels as functions of radius
	if(ider.ne.0) then
c------------------------------LAPO(10/98)------------------------------
	DQDVPH(jpath,i12)=qints(5)-sumdqdvph
	sumdqdvph=qints(5)
	DQDVPv(jpath,i12)=qints(9)-sumdqdvpv
	sumdqdvpv=qints(9)
	DQDVsH(jpath,i12)=qints(13)-sumdqdvsh
	sumdqdvsh=qints(13)
	DQDVsv(jpath,i12)=qints(17)-sumdqdvsv
	sumdqdvsv=qints(17)
	DQDeta(jpath,i12)=qints(21)-sumdqdeta
	sumdqdeta=qints(21)

c	onewayker1(jpath)=aktemp(1)
c	onewayker2(jpath)=aktemp(2)
c	onewayker3(jpath)=aktemp(3)
c	onewayker4(jpath)=aktemp(4)
c	onewayker5(jpath)=aktemp(5)
c	rayker(jpath)=0.5d0*(x1+x2)+xnorm
c------------------------------LAPO(10/98)------------------------------
	endif
c---------------------------------LAPO----------------------------------
c---------------------------------LAPO----------------------------------
c	write(124,*)jpath,n,((qints(1)-sumd)*180.)/pi,
c     &x1+xnorm,x2+xnorm
c---------------------------------LAPO----------------------------------

        sumd=qints(1)

c---------------------------------LAPO----------------------------------
c--save incremental travel time as a function of radius
	timearray(jpath,i12)=qints(3)-sumti
	sumti=qints(3)
c---------------------------------LAPO----------------------------------

c--the next is a SIGNIFICANT change that I made (8-26-98) 
c---------------------------------old-----------------------------------
c        if(x1+xnorm.gt.rdep) kdep=jpath
c---------------------------------old-----------------------------------
c---------------------------------LAPO----------------------------------
        if(x2+xnorm.eq.rdep)then
	   kdep=jpath
c	   write(124,*)'kdep=',jpath
	endif
c---------------------------------LAPO----------------------------------

        goto 10
   20   continue
c ----------------------------------------------------------------------
      if(itim.eq.1) return
      if(xnorm.eq.0.d0) return
      if(ifl.ge.0) qints(2)=qints(2)-dxtdp*fxtpp/(dsqrt(xt(n)-xnorm))
      return
      end

      subroutine qtau(x,iq,q1)
c  returns the radicand of qtau for radius x
      implicit double precision (a-h,o-z)
      logical isotrp
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/parm/rho,vpv,vsv,vph,vsh,eta,s(5),r
      common/parm1/xa,xc,xn,xl,xf,tdif(5)
      common/isot/isotrp
      y=x/rnorm  
      call getmod(y,iq,rho,vpv,vsv,vph,vsh,eta)
      xc=rho*vpv*vpv
      xl=rho*vsv*vsv
      xa=xc
      if(.not.isotrp) then
        xa=rho*vph*vph
        xn=rho*vsh*vsh
        xf=eta*(xa-2.d0*xl)
        if(imod.ne.0) then
          xkapa=(4.d0*xa+xc+4.d0*xf-4.d0*xn)/9.d0
          xmu=(xa+xc-2.d0*xf+5.d0*xn+6.d0*xl)/15.d0
          xa=xkapa+4.d0*xmu/3.d0
          xc=xa  
          xf=xkapa-2.d0*xmu/3.d0
          xn=xmu 
          xl=xn
          vph=dsqrt(xa/rho)
          vpv=vph
          vsh=dsqrt(xn/rho)
          vsv=vsh
          eta=1.d0
        end if
      end if
c*** avoid the center of the earth!
      if(x.lt.1.d-3) x=1.d-3
      px2=(p/x)**2
      if(itype.eq.1) then
c**** SH
c  reflect at fluid boundary
        if(vsv.eq.0.d0) then
          q1=-1.d0
          return
        end if
        if(isotrp) then
          q1=1.d0/(vsv*vsv)-px2
        else 
          q1=(1.d0-vsh*vsh*px2)/(vsv*vsv)
        end if
      else if(vsv.eq.0.d0) then
c  reflect SV at fluid boundary
        if(itype.eq.3) then
           q1=-1.d0
           return
        end if
c*** P in fluid
        q1=1.d0/(vpv*vpv)-px2
      else
c*** P/SV in solid
        if(isotrp) then
          if(itype.eq.2) q1=1.d0/(vpv*vpv)-px2
          if(itype.eq.3) q1=1.d0/(vsv*vsv)-px2
        else
          aa=0.5d0/(vsv*vsv)
          bb=0.5d0/(vpv*vpv)
          s(1)=aa+bb
          s(2)=aa-bb
          s(3)=(xa*xc-xf*xf-2.d0*xf*xl)/(2.d0*xc*xl)
          s(4)=s(3)*s(3)-xa/xc
          s(5)=0.5d0*(1.d0+xa/xl)/(vpv*vpv)-s(1)*s(3)
          r=sqrt((s(4)*px2+2.d0*s(5))*px2+s(2)*s(2))
          if(itype.eq.2) q1=s(1)-s(3)*px2-r
          if(itype.eq.3) q1=s(1)-s(3)*px2+r
        end if
      end if
      return
      end

      subroutine dxt(x,n,d,f)
c     accepts as input ray parameter (through common/in/),
c     turning radius x (real or fictitious), and layer
c     number n of layer containing x.
c     returns the derivative d of the turning radius w.r.t. p
c
c     the f.p. in the calling routine (integr) are:
c              x -- xnorm
c              n -- n (layer number)
c              d -- dxtdp
c              f -- fxttp
c
      implicit double precision (a-h,o-z)
      logical isotrp
      dimension v(3),dv(3),ds(5),bigv(3,3),bigd(3,3)
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/parm/rho,vpv,vsv,vph,vsh,eta,s(5),r
      common/parm1/xa,xc,xn,xl,xf,tdif(5)
      common/isot/isotrp
c
      y=x/rnorm
      call getmod(y,n,rho,vpv,vsv,vph,vsh,eta)
      xc=rho*vpv*vpv
      xl=rho*vsv*vsv
      if(.not.isotrp) then
        xa=rho*vph*vph
        xn=rho*vsh*vsh
        xf=eta*(xa-2.d0*xl)
        if(imod.ne.0) then
          xkapa=(4.d0*xa+xc+4.d0*xf-4.d0*xn)/9.d0
          xmu=(xa+xc-2.d0*xf+5.d0*xn+6.d0*xl)/15.d0
          xa=xkapa+4.d0*xmu/3.d0
          xc=xa
          xf=xkapa-2.d0*xmu/3.d0
          xn=xmu
          xl=xn
          vph=dsqrt(xa/rho)
          vpv=dsqrt(xc/rho)
          vsh=dsqrt(xn/rho)
          vsv=dsqrt(xl/rho)
          eta=1.d0
        end if
      end if
      d=0.
      f=0.
      if(x.lt.1.d-3) return
      px2=(p/x)**2
c
      if(itype.eq.1) then
        if(vsv.eq.0.d0) return
c**** SH
        if(.not.isotrp) then
          dvsh=cder(y,7,n)/vsh
          dvsv=cder(y,3,n)/vsv
          top=px2*xn/(xl*p)
          bot=(px2/x-dvsv/(vsh*vsh)+px2*(dvsv-dvsh))*xn/xl
        else
          top=px2/p
          bot=px2/x-cder(y,3,n)/(vsv**3)
        end if
        d=top/bot
        f=0.5d0*d*dsqrt(2.d0*bot)
      else if(vsv.eq.0.d0) then
c**** fluid
        top=px2/p
        bot=px2/x-cder(y,2,n)/(vpv**3)
        d=top/bot
        f=0.5d0*d*dsqrt(2.d0*bot)
      else
c**** P/SV
        v(2)=vpv
        v(3)=vsv
        dv(2)=cder(y,2,n)
        dv(3)=cder(y,3,n)
        a=dv(2)/(v(2)**3)
        b=dv(3)/(v(3)**3)
        ds(1)=-a-b
        ds(2)=a-b
        if(isotrp) then
          top=2.d0*px2/p
          bot=ds(1)+2.d0*px2/x
          if(itype.eq.2) bot=bot-ds(2)
          if(itype.eq.3) bot=bot+ds(2)
        else
          v(1)=vph
          dv(1)=cder(y,6,n)
          deta=cder(y,8,n)
          do 1 i=1,3
          do 1 j=1,3
    1       bigv(i,j)=v(i)*v(i)/(v(j)*v(j))
          do 2 i=1,3
          do 2 j=1,3
    2       bigd(i,j)=2.d0*bigv(i,j)*(dv(i)/v(i)-dv(j)/v(j))
          aa=0.5d0/(vsv*vsv)
          bb=0.5d0/(vpv*vpv)
          s(1)=aa+bb
          s(2)=aa-bb
          s(3)=(xa*xc-xf*xf-2.d0*xf*xl)/(2.d0*xc*xl)
          s(4)=s(3)*s(3)-xa/xc
          s(5)=0.5d0*(1.d0+xa/xl)/(vpv*vpv)-s(1)*s(3)
          r=sqrt((s(4)*px2+2.d0*s(5))*px2+s(2)*s(2))
c*** ds(1--5) are derivs of s() wrt x (similarly with rdx and r)
          b1=-0.5d0*eta*bigd(1,3)+deta*(2.d0-0.5d0*bigv(1,3))
          a2=-1.d0+eta*(2.d0-0.5d0*bigv(1,3))
          b2=eta*bigd(1,2)+deta*bigv(1,2)
          a3=-2.d0*eta*deta*bigv(3,2)
          a4=2.d0*(1.d0-eta)*(eta*bigd(3,2)+deta*bigv(3,2))
          ds(3)=0.5d0*bigd(1,3)+eta*bigv(1,2)*b1+a2*b2+a3+a4
          ds(4)=2.d0*s(3)*ds(3)-bigd(1,2)
          ds(5)=(0.5d0*bigd(1,3)-dv(2)*(1.d0+bigv(1,3))/vpv)/(vpv*vpv)
     +          -s(1)*ds(3)-s(3)*ds(1)
          a=0.5d0*px2*px2*(ds(4)-4.d0*s(4)/x)
          b=px2*(ds(5)-2.d0*s(5)/x)+s(2)*ds(2)
          rdx=(a+b)/r
          rdp=2.d0*px2*(s(4)*px2+s(5))/(r*p)
          top=2.d0*s(3)*px2/p
          bot=ds(1)+px2*(2.d0*s(3)/x-ds(3))
          if(itype.eq.2)then
            top=top+rdp
            bot=bot-rdx
          else
            bot=bot+rdx
            top=top-rdp
          end if
        end if
        d=top/bot
        f=0.5d0*d*dsqrt(bot)
      end if
      return
      end

      subroutine fqanis(x,iq,vals)
c   builds up the discretized integrands (with appropriate
c   renormalization) for use in subroutine integr.
c              x = modified radius
c             iq = layer number
c         vals() = returned values of desired integrands
c
c  (rewritten to handle vertical incidence, p=0)
      implicit double precision (a-h,o-z)
      dimension vals(*)
c	DIMENSION zk(4)

	dimension refe(5)

      logical isotrp
      common/isot/isotrp
      common/parm/rho,vpv,vsv,vph,vsh,eta,s(5),r
      common/parm1/xa,xc,xn,xl,xf,tdif(5)
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/xnorm/xnorm
      common/dx/dxtdp,fxtpp,ifl
c*** here, y is true radius, xnorm is turning point and x=sqrt(y-xnorm) 
c*** so y=xnorm+x*x and dy=2xdx, this change of variable eliminates 
c*** turning point singularities
      y=x*x+xnorm
      call qtau(y,iq,q)
      q=dsqrt(q)
      p2=p/(y*y)
      py2=p2*p
c*** delta
      q1=p2/q
      f=1.d0
      if(.not.isotrp) then
        if(itype.eq.1) then
          f=xn/xl
        else
          if(vsv.ne.0.d0) then
            a=(s(4)*py2+s(5))/r
            if(itype.eq.3) a=-a
            f=s(3)+a
          end if
        end if
      end if
      q1=q1*f
      vals(1)=2.d0*x*q1
      if(itim.eq.1) return
c*** dX/dp
      q2=f*(1.d0+p*q1/q)/(y*y*q)
      if(.not.isotrp) then
        if(itype.ne.1)then
          a=2.d0*py2*(s(4)*s(2)*s(2)-s(5)*s(5))/(q*(r**3)*(y**2))
          if(itype.eq.3) a=-a
          q2=q2+a
        end if
      end if
      if(ifl.ge.0) q2=q2-fxtpp*dxtdp/(2.d0*x**3)
      vals(2)=2.d0*x*q2
      if(itim.eq.2) return
c*** time and tstar     
      q3=q+p*q1
      vals(3)=2.d0*x*q3
      z=y/rnorm
      call getq(z,iq,qalpha,qbeta)
      if(itype.eq.2) then
        vals(4)=vals(3)*qalpha
      else
        vals(4)=vals(3)*qbeta
      end if
      if(ider.eq.0) return
c*** derivs
c  the tdifs are the derivs wrt vph,vpv,vsh,vsv,eta
      call tder(y,q)
c  we want to integrate these over each layer multiplied
c  by the appropriate powers of radius for a prem 
c  parameterization -- this will need to be changed for
c  other parametizations- YES! (10/98)
C--At this stage I just compute the increment in the radial integral.
C--in analogy with what is done with travel time.
c-------------------------------old(10/98)------------------------------
c      zk(1)=1.d0
c      do k=2,4
c        zk(k)=z*zk(k-1)
c      enddo
c-------------------------------old(10/98)------------------------------
c------------------------------LAPO(10/98)------------------------------
c--the prem-velocities here are in km/s.
	refe(1)=vph
	refe(2)=vpv
	refe(3)=vsh
	refe(4)=vsv
	refe(5)=eta
c	write(46,*)6371.-y,vpv
c	write(47,*)6371.-y,vsv
c	write(48,*)6371.-y,vph
c	write(49,*)6371.-y,vsh
c------------------------------LAPO(10/98)------------------------------

      do 10 j=1,5

c------------------------------LAPO(10/98)------------------------------
	val0=refe(j)
c------------------------------LAPO(10/98)------------------------------

      do 10 k=1,4
        ind=4*(j-1)+k+4

c------------------------------LAPO(10/98)------------------------------
c   10   vals(ind)=2.d0*x*tdif(j)*zk(k)
c--I also need to multiply the kernels by the PREM reference values 
c--val0, so that the unknown of the inverse problem is defined as a 
c--RELATIVE perturbation.
   10   vals(ind)=2.d0*x*tdif(j)*val0
c------------------------------LAPO(10/98)------------------------------

      return
      end

      subroutine tder(x,qq)
c*** computes the derivatives of qtau wrt vph,vpv,vsh,vsv,eta
c*** at constant p, they are returned in array tdif in common parm1
c*** note that qtau must have been called first and qq=qtau is input
      implicit double precision (a-h,o-z)
      logical isotrp
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/parm/rho,vpv,vsv,vph,vsh,eta,s(5),r
      common/parm1/xa,xc,xn,xl,xf,tdif(5)
      common/isot/isotrp
      dimension sd(5,5)
      do i=1,5
        tdif(i)=0.d0
      enddo
      if(x.lt.1.d-3) return
      px2=(p/x)**2
      if(itype.eq.1) then
c**** SH 
        if(isotrp) then
          tdif(4)=-1.d0/(qq*vsv**3)
          tdif(3)=tdif(4)
        else
          tdif(3)=-px2*vsh/(qq*vsv*vsv)
          tdif(4)=-qq/vsv
        end if
      else if(vsv.eq.0.d0) then
c*** P in fluid
        tdif(2)=-1.d0/(qq*vpv**3)
        tdif(1)=tdif(2)
      else
c*** P/SV in solid
        if(isotrp) then
          if(itype.eq.2) then
            tdif(1)=-1.d0/(qq*vpv**3)
            tdif(2)=tdif(1)
          else
            tdif(3)=-1.d0/(qq*vsv**3)
            tdif(4)=tdif(3)
          end if
        else
c*** note that sd() are derivs of s() wrt A,C,N,L,F resp
c*** similarly rd is deriv of r
          sd(1,1)=0.d0
          sd(1,2)=-rho/(2.d0*xc*xc)
          sd(1,4)=-rho/(2.d0*xl*xl)
          sd(1,5)=0.d0
          sd(2,1)=0.d0
          sd(2,2)=-sd(1,2)
          sd(2,4)=sd(1,4)
          sd(2,5)=0.d0
          sd(3,1)=1.d0/(2.d0*xl)
          sd(3,2)=(xf*xf+2.d0*xl*xf)/(2.d0*xl*xc*xc)
          sd(3,4)=(xf*xf-xa*xc)/(2.d0*xc*xl*xl)
          sd(3,5)=-(1.d0+xf/xl)/xc
          sd(4,1)=s(3)/xl-1.d0/xc
          sd(4,2)=2.d0*s(3)*sd(3,2)+xa/(xc*xc)
          sd(4,4)=2.d0*s(3)*sd(3,4)
          sd(4,5)=2.d0*s(3)*sd(3,5)
          sd(5,1)=rho/(2.d0*xl*xc)-s(1)*sd(3,1)
          sd(5,2)=-rho*(1.d0+xa/xl)/(2.d0*xc*xc)
          sd(5,2)=sd(5,2)-s(1)*sd(3,2)-s(3)*sd(1,2)
          sd(5,4)=-rho*xa/(2.d0*xc*xl*xl)-s(1)*sd(3,4)-s(3)*sd(1,4)
          sd(5,5)=-s(1)*sd(3,5)
          do j=1,5
            if(j.ne.3) then
              rd=((px2*sd(4,j)+2.d0*sd(5,j))*px2+
     +                2.d0*s(2)*sd(2,j))/(2.d0*r)
              if (itype.eq.3) rd=-rd
              tdif(j)=(sd(1,j)-sd(3,j)*px2-rd)/(2.d0*qq)
            end if
          enddo

c--I checked the following formulae and I think they are ok...
c*** now we change from A,C,N,L,F to A,C,N,L,eta
          tdif(1)=tdif(1)+eta*tdif(5)
          tdif(4)=tdif(4)-2.d0*eta*tdif(5)
          tdif(5)=tdif(5)*(xa-2.d0*xl)

c*** now we change to vph,vpv,vsh,vsv,eta
          tdif(1)=2.d0*rho*vph*tdif(1)
          tdif(2)=2.d0*rho*vpv*tdif(2)
          tdif(4)=2.d0*rho*vsv*tdif(4)

        end if
      end if
      return
      end

      double precision function f(pp)
c   f is delta as a function of p if itim=1 in calling routine
c   f is d(delta)/d(p) as a function of p if itim=2
      implicit double precision (a-h,o-z)
      dimension qvec(500)
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      call deriv(pp,qvec,xtu)
      f=qvec(itim)
      return
      end

      subroutine gauslv(r1,r2,nin,f,iq,fint,nint)
c action: performs gauss-legendre integrations
c input parameters (supplied by routine integr):
c     r1,r2  - integration bounds
c     nin    - (idleg in calling routine) order of the legendre polynomial
c     f(external)- fqanis (integral kernel)
c     iq     - number of the layer
c     fint   - array, contains the values (delta,t,d(delta)/d(p) and
c              partials derivatives) requested by the user (qints(3)=
c              fint(3) will undergo a later scaling in routine integr)
c     nint   - number of elements of fint requested by the user
c  ***** nint must be less than nmax *****
c
c the choice of the order of the legendre integration is limited to
c a few values (recorded via a data statement in array nord). through
c nnord and nord the nearest higher selected order is choosen. maximum
c order is 24. array x contains the zeros and w contains the weights.
c array iad contains the addresses of the first zero and weight for
c each polynomial order. the symmetry of zeros and weights is welcomed.
c
c    the integration bounds are r1,r2. they are transformed to -1,1
c with the new variable t: x=((r1+r2)/2)+t*((r2-r1)/2). the integral
c is computed in t and then multiplied by (r2-r1)/2.
c
      implicit double precision(a-h,o-z)
      parameter (nmax=500)
      dimension fint(nint),vals(nmax),vals1(nmax),sum(nmax)
      dimension iad(13),nord(24),x(65),w(65),w1(29),w2(36)
     1     ,x1(29),x2(36),nnord(13)
      equivalence (x(1),x1(1)),(x(30),x2(1)),(w(1),w1(1))
     1     ,(w(30),w2(1))
      data  iad/ 1, 2, 4, 6, 9,12,16,20,25,30,36,44,54/
      data nord/1,1,2,3,4,5,6,7,8,9,10,10,11,11,11,11
     1     ,12,12,12,12,13,13,13,13/
      data nnord/ 2, 3, 4, 5, 6, 7, 8, 9,10,12,16,20,24/
      data x1/ .57735 02691 89626d0,.00000 00000 00000d0
     1        ,.77459 66692 41483d0,.33998 10435 84856d0
     2        ,.86113 63115 94053d0,.00000 00000 00000d0
     3        ,.53846 93101 05683d0,.90617 98459 38664d0
     4        ,.23861 91860 83197d0,.66120 93864 66265d0
     5        ,.93246 95142 03152d0,.00000 00000 00000d0
     6        ,.40584 51513 77397d0,.74153 11855 99394d0
     7        ,.94910 79123 42759d0,.18343 46424 95650d0
     8        ,.52553 24099 16329d0,.79666 64774 13627d0
     9        ,.96028 98564 97536d0,.00000 00000 00000d0
     x        ,.32425 34234 03809d0
     1        ,.61337 14327 00590d0,.83603 11073 26636d0
     2        ,.96816 02395 07626d0,.14887 43389 81631d0
     3        ,.43339 53941 29247d0,.67940 95682 99024d0
     4        ,.86506 33666 88985d0,.97390 65285 17172d0/
      data x2/ .12523 34085 11469d0,.36783 14989 98180d0
     1        ,.58731 79542 86617d0,.76990 26741 94305d0
     2        ,.90411 72563 70475d0,.98156 06342 46719d0
     3        ,.09501 25098 37637d0,.28160 35507 79259d0
     4        ,.45801 67776 57227d0,.61787 62444 02644d0
     5        ,.75540 44083 55003d0,.86563 12023 87831d0
     6        ,.94457 50230 73233d0,.98940 09349 91650d0
     7        ,.07652 65211 33497d0,.22778 58511 41645d0
     8        ,.37370 60887 15420d0,.51086 70019 50827d0
     9        ,.63605 36807 26515d0,.74633 19064 60151d0
     x        ,.83911 69718 22218d0,.91223 44282 51325d0
     1        ,.96397 19272 77913d0,.99312 85991 85094d0
     2        ,.06405 68928 62605d0,.19111 88674 73616d0
     3        ,.31504 26796 96163d0,.43379 35076 26045d0
     4        ,.54542 14713 88840d0,.64809 36519 36976d0
     5        ,.74012 41915 78554d0,.82000 19859 73902d0
     6        ,.88641 55270 04401d0,.93827 45520 02733d0
     7        ,.97472 85559 71309d0,.99518 72199 97021d0/
      data w1/1.00000 00000 00000d0,.88888 88888 88889d0
     1        ,.55555 55555 55556d0,.65214 51548 62546d0
     2        ,.34785 48451 37454d0,.56888 88888 88889d0
     3        ,.47862 86704 99366d0,.23692 68850 56189d0
     4        ,.46791 39345 72691d0,.36076 15730 48139d0
     5        ,.17132 44923 79170d0,.41795 91836 73469d0
     6        ,.38183 00505 05119d0,.27970 53914 89277d0
     7        ,.12948 49661 68870d0,.36268 37833 78362d0
     8        ,.31370 66458 77887d0,.22238 10344 53374d0
     9        ,.10122 85362 90376d0,.33023 93550 01260d0
     x        ,.31234 70770 40003d0
     1        ,.26061 06964 02935d0,.18064 81606 94857d0
     2        ,.08127 43883 61574d0,.29552 42247 14753d0
     3        ,.26926 67193 09996d0,.21908 63625 15982d0
     4        ,.14945 13491 50581d0,.06667 13443 08688d0/
      data w2/ .24914 70458 13403d0,.23349 25365 38355d0
     1        ,.20316 74267 23066d0,.16007 83285 43346d0
     2        ,.10693 93259 95318d0,.04717 53363 86512d0
     3        ,.18945 06104 55068d0,.18260 34150 44923d0
     4        ,.16915 65193 95003d0,.14959 59888 16576d0
     5        ,.12462 89712 55534d0,.09515 85116 82493d0
     6        ,.06225 35239 38648d0,.02715 24594 11754d0
     7        ,.15275 33871 30726d0,.14917 29864 72604d0
     8        ,.14209 61093 18382d0,.13168 86384 49176d0
     9        ,.11819 45319 61518d0,.10193 01198 17240d0
     x        ,.08327 67415 76705d0,.06267 20483 34109d0
     1        ,.04060 14298 00387d0,.01761 40071 39152d0
     2        ,.12793 81953 46752d0,.12583 74563 46828d0
     3        ,.12167 04729 27803d0,.11550 56680 53726d0
     4        ,.10744 42701 15966d0,.09761 86521 04114d0
     5        ,.08619 01615 31953d0,.07334 64814 11080d0
     6        ,.05929 85849 15437d0,.04427 74388 17420d0
     7        ,.02853 13886 28934d0,.01234 12297 99987d0/
c
      if(nint.gt.nmax) then
         print *,'nmax in routine gauslv needs increasing',nint,nmax
         stop
      end if
      n=nin
      if(nin.lt.1) n=1
      if(nin.gt.24) n=24
c ----nin is the requested degree of the integration
c ----the actual order is nnord(nord(n)); the zeros are from
c ----x(iad(nord(n)) on, the weights start at w(iad(nord(n)))=w(ia)
      ind=nord(n)
      n=nnord(ind)
      ia=iad(ind)
      nc=n/2
      nc2=2*nc
      y1=.5d0*(r2+r1)
      y2=.5d0*(r2-r1)
      do j=1,nint
        sum(j)=0.d0
      enddo
      if(n.ne.nc2) then
        call f(y1,iq,vals)
        do j=1,nint
          sum(j)=w(ia)*vals(j)
        enddo
        ia=1+ia
      end if
      if(nc.ne.0) then
        do 10 i=1,nc
        t1=x(ia)*y2
        call f(y1+t1,iq,vals)
        call f(y1-t1,iq,vals1)
        do j=1,nint
          sum(j)=sum(j)+w(ia)*(vals(j)+vals1(j))
        enddo
   10   ia=1+ia
      end if
      do j=1,nint
        fint(j)=fint(j)+y2*sum(j)
      enddo
      return
      end

      subroutine getmod(y,iq,rho,vpv,vsv,vph,vsh,eta)
c*** returns model vector at normalized radius y which
c*** shell iq. When iq is changed, the routine checks
c*** to see if the layer is isotropic (returned through 
c*** common isot). 
      implicit real*8(a-h,o-z)
      logical isotrp
      dimension test(4)
      common/coeff/coef(4,8,20)
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/isot/isotrp
      data iqs/-1/,test/1.d0,0.d0,0.d0,0.d0/
      if(iq.ne.iqs) then
        iqs=iq
        isotrp=.false.
        if(iso.ne.1) then
          sum=0.d0
          do i=1,4
            sum=sum+coef(i,2,iq)-coef(i,6,iq)
     +        +coef(i,3,iq)-coef(i,7,iq)+coef(i,8,iq)-test(i)
          enddo
          if(sum.eq.0.d0) isotrp=.true.
        end if
      end if
      rho=coef(1,1,iq)+y*(coef(2,1,iq)+y*(coef(3,1,iq)+y*coef(4,1,iq)))
      vpv=coef(1,2,iq)+y*(coef(2,2,iq)+y*(coef(3,2,iq)+y*coef(4,2,iq)))
      vsv=coef(1,3,iq)+y*(coef(2,3,iq)+y*(coef(3,3,iq)+y*coef(4,3,iq)))
      if(.not.isotrp) then
      vph=coef(1,6,iq)+y*(coef(2,6,iq)+y*(coef(3,6,iq)+y*coef(4,6,iq)))
      vsh=coef(1,7,iq)+y*(coef(2,7,iq)+y*(coef(3,7,iq)+y*coef(4,7,iq)))
      eta=coef(1,8,iq)+y*(coef(2,8,iq)+y*(coef(3,8,iq)+y*coef(4,8,iq)))
      else
        vph=vpv
        vsh=vsv
        eta=1.d0
      end if
      return
      end

      subroutine getq(y,iq,qalpha,qbeta)
c*** returns qalpha and qbeta at normalized radius y which
c*** shell iq. Note that this assumes that the 5 and 4
c*** model coefficients relate to inverse Qmu and inverse Qkappa
c*** which is the case for aniprmc
c*** this code is a little cavalier about anisotropic layers!
c*** also we assume constant Q in layers (as is done in readmod)
      implicit real*8(a-h,o-z)
      common/coeff/coef(4,8,20)
      vpv=coef(1,2,iq)+y*(coef(2,2,iq)+y*(coef(3,2,iq)+y*coef(4,2,iq)))
      vsv=coef(1,3,iq)+y*(coef(2,3,iq)+y*(coef(3,3,iq)+y*coef(4,3,iq)))
      qka=coef(1,4,iq)
      qmu=coef(1,5,iq)
      qbeta=qmu
      qalpha=qka+4.d0/3.d0*(qmu-qka)*(vsv/vpv)**2
      return
      end

      double precision function cder(y,m,iq)
c*** finds radial derivative of model parameter 'm' where
c  m=1 rho; m=2 vpv, =3 vsv; =4 Qm; =5 Qk; =6 vph, =7 vsh; =8 eta
      implicit real*8(a-h,o-z)
      common/coeff/coef(4,8,20)
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      cder=(coef(2,m,iq)+y*(2.d0*coef(3,m,iq)+3.d0*coef(4,m,iq)
     +      *y))/rnorm
      return
      end

      subroutine readmod(filnam)
c     read in the model
c         coef(4,8,20) contains the 4 polynomial coefficients for
c                 each of 8 parameters for up to 20 layers
c*** assume that first layer is ic and second layer is oc
      implicit real*8(a-h,o-z)
      character*256 filnam,modnam
      common/coeff/coef(4,8,20)
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/slow/picblp,picbls,picbup,pcmblp,pcmbup,pcmbus,
     +   p660lp,p660ls,lic,loc,llm,lsrc,klay(20)

c------------------------------Lapo(11/98)------------------------------
	common/SLOWLAPO/pcmblpLAPO,pcmbupLAPO
c------------------------------Lapo(11/98)------------------------------

      data pi/3.14159265358979d0/,rtol/0.d0/
      rad=180.d0/pi

c      print *,'enter model file name :'
c      read(*,'(a256)') filnam

      open(8,file=filnam)
      read(8,'(a256)') modnam
      read(8,*) ifanis,tref,ifdeck
      read(8,*) numlyr,nic,noc,rnorm
c      read(8,900) numlyr,nic,noc,rnorm,moho
c      write(99,800) numlyr,nic,noc,rnorm,moho
c  800 format(3i5,f10.3,i5)
c  900 format(3i5,d15.9,i5)
      npar=5
      if(ifanis.ne.0) npar=8
      do 5 i=1,numlyr
        k=numlyr-i+1
        read(8,*) nl(k),xb(k),xt(k)
c        read(8,905) nl(k),junk,xb(k),xt(k)
c        write(99,805) nl(k),xb(k),xt(k)
c  805 format(i5,2f10.3)
c  905   format(2i5,2d15.9)
c        dum=0.d0
        do j=1,npar
          read(8,810)(coef(jj,j,k),jj=1,4)
c*** this assumes polynomials are in big Q and are constant only
          if(j.eq.4.or.j.eq.5.and.coef(1,j,k).ne.0.d0)
     +                          coef(1,j,k)=1.d0/coef(1,j,k)
c          write(99,820)(coef(jj,j,k),jj=1,4),dum
  810 format(5g16.9)
c  820 format(5g16.9)
c  810 format(5g9.5)
c  910   format(4d16.9)
        enddo
c*** i'm not sure if i really need to do this
        if(npar.eq.5) then
          do kk=1,4
            coef(kk,6,k)=coef(kk,2,k)
            coef(kk,7,k)=coef(kk,3,k)
            coef(kk,8,k)=0.d0
          enddo
          coef(1,8,k)=1.d0
        end if
        if(abs(xt(k)-5700.).lt.20.) llm=k
    5   continue
c**** we are assuming a prem-type parameterization for the derivs
c  ie, cubic polynomials so k=1,4 and we have 5 parameters (vph,
c  vpv,vsh,vsv,eta) so we have 20 coefficients per layer
      nplay=20
      close(8)
c  find slownesses just below and above icb,cmb, and 660      
      lic=numlyr
      loc=numlyr-1
      iq=numlyr
      y=(xt(iq)-rtol)/rnorm
      call getmod(y,iq,rho,vpv,vsv,vph,vsh,eta)
      picblp=y*rnorm/vpv/rad
      picbls=y*rnorm/vsv/rad
      iq=numlyr-1
      y=(xb(iq)+rtol)/rnorm
      call getmod(y,iq,rho,vpv,vsv,vph,vsh,eta)
      picbup=y*rnorm/vpv/rad
      y=(xt(iq)-rtol)/rnorm
      call getmod(y,iq,rho,vpv,vsv,vph,vsh,eta)
      pcmblp=y*rnorm/vpv/rad

c------------------------------Lapo(11/98)------------------------------
	pcmblpLAPO=1.d0/vpv
c------------------------------Lapo(11/98)------------------------------

      iq=numlyr-2
      y=(xb(iq)+rtol)/rnorm
      call getmod(y,iq,rho,vpv,vsv,vph,vsh,eta)
      pcmbup=y*rnorm/vpv/rad

c------------------------------Lapo(11/98)------------------------------
	pcmbupLAPO=1.d0/vpv
c------------------------------Lapo(11/98)------------------------------

      pcmbus=y*rnorm/vsv/rad
      iq=llm
      y=(xt(iq)-rtol)/rnorm
      call getmod(y,iq,rho,vpv,vsv,vph,vsh,eta)
      p660lp=y*rnorm/vpv/rad
      p660ls=y*rnorm/vsv/rad
      return
      end

      subroutine raynam(ray,ierr)
      implicit real*8(a-h,o-z)
C only names using P,S,I,J,K,p,s,i,c allowed
c p, P, s, or S must be first character.
c The reflection and turning ray , e.g. P and PcP are distinguished 
c by the ray parameter range
c
      common/ray$/pmin,pmax,nps(20,2),npsrc(2),nn,iup,jref,npleg,nsleg
      common/slow/picblp,picbls,picbup,pcmblp,pcmbup,pcmbus,
     +   p660lp,p660ls,lic,loc,llm,lsrc,klay(20)
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      logical split
      character*(*) ray      
      data ptol/1.d-6/
c*** split is false if source is at existing interface (eg surface)
c*** if split is true, we put multiplicity for legs above source in
c*** npsrc
      split=.true.
      if(rdep.eq.xt(lsrc)) split=.false.
c*** get indices for bottom of layers -- note
c*** that the source layer and the turning layer will 
c*** have two and one  extra points resp.
      klay(1)=nl(1)-1

c------------------------SIGNIFICANT CHANGE (8-27-98)
C      if(lsrc.eq.1.and.split) klay(1)=klay(1)+2
      if(lsrc.eq.1.and.split) klay(1)=klay(1)*2
c----------------------------------------------------

      do i=2,numlyr
        klay(i)=nl(i)+klay(i-1)-1

c------------------------SIGNIFICANT CHANGE (8-27-98)
C        if(lsrc.eq.i.and.split) klay(i)=klay(i)+2
        if(lsrc.eq.i.and.split) klay(i)=(nl(i)*2)+klay(i-1)-2
c----------------------------------------------------

      enddo
      ierr=0                                
c*** iup=1 if little p or 2 if little s
      iup=0
c*** jref=1 means reflection from cmb, 2= reflection from icb
      jref=0
      npsrc(1)=0
      npsrc(2)=0
      do 5 i=1,numlyr
        nps(i,1)=0
    5   nps(i,2)=0                                                    
C ************************************************* evaluate # ray legs
      j1=0
      do 10 i=1,70
      j=ichdec(ray(i:i),k)
      goto(10,11,11,12,12,13,14,14,15,16,20),j
C**** P or S ****
   11 goto(20,111,112,20,20,114,113,113,20,100),j1
      call addray(nps(1,k),lsrc,loc-1,1)
      goto 100
  111 call addray(nps(1,1),1,loc-1,1)              
      if(split) npsrc(1)=npsrc(1)+1
      call addray(nps(1,k),1,loc-1,1)
      if(split) npsrc(k)=npsrc(k)+1
      goto 100
  112 call addray(nps(1,2),1,loc-1,1)              
      if(split) npsrc(2)=npsrc(2)+1
  113 call addray(nps(1,k),1,loc-1,1)              
      if(split) npsrc(k)=npsrc(k)+1
      goto 100
  114 call addray(nps(1,1),loc,lic-1,1)            
      goto 100
C**** I or J ****
   12 goto(20,20,20,121,121,121,20,20,20,20),j1
  121 call addray(nps(1,k),lic,numlyr,2)           
      goto 100
C**** K ****
   13 goto(20,132,132,100,100,131,20,20,100,20),j1
  131 call addray(nps(1,1),loc,lic-1,1)            
  132 call addray(nps(1,1),loc,lic-1,1)            
      goto 100
C**** p or s ****
   14 if(j1.ne.0) goto 20
      iup=k
      call addray(nps(1,k),1,lsrc-1,1)
      if(split) npsrc(k)=npsrc(k)+1
      goto 100
c**** i ****
   15 jref=2
      if(j1.ne.6) goto 20
      goto 100
C**** c ****
   16 jref=1
      if(j1.lt.2.or.j1.gt.3) goto 20
  100 j1=j  
      k1=k
   10 continue
C****  receiver leg ****
      goto(20,201,201,20,20,20,201,201,20,20),j1
  201 call addray (nps(1,k1),1,loc-1,1)
      if(split) npsrc(k1)=npsrc(k1)+1
c***************************************************** finish of ray leg summation
c*** nn is total # layers sampled
      do j=numlyr,1,-1
        nn=j
        if(nps(j,1)+nps(j,2).ne.0) goto 300
      enddo
  300 npleg=0
      nsleg=0
      do i=1,nn
        npleg=npleg+nps(i,1)
        nsleg=nsleg+nps(i,2)
      enddo
c*** evaluate ray parameter range
      if(jref.eq.1.and.nn.lt.loc) then
c*** bounce off core (no legs beneath unlike e.g. PcPPKP)
        pmax1=1.e10
        pmax2=1.e10
        if(nps(loc-1,1).ne.0) pmax1=pcmbup
        if(nps(loc-1,2).ne.0) pmax2=pcmbus
        pmax=min(pmax1,pmax2)-ptol
        pmin=0.
c*** bounce of inner core -- no legs beneath
      else if(jref.eq.2.and.nn.lt.lic) then
        pmax=picbup-ptol
        pmin=0.
      else
c*** turning rays
c*** following for rays turning in outer core
        if(nn.eq.loc) then
          pmin=picbup+ptol
          pmax1=1.e10
          pmax2=1.e10
          if(nps(loc-1,2).ne.0) pmax1=pcmbus
          if(nps(loc-1,1).ne.0) pmax1=pcmbup
          pmax=min(pmax1,pmax2)
          pmax=min(pmax,pcmblp)-ptol
          return
c*** following for rays turning in inner core
c*** program has numerical problems if we turn very close
c  to center of Earth ... hence limit on pmin
        else if(nn.eq.lic) then
          pmin=.0001
          pmax1=1.e10
          pmax2=1.e10
          if(nps(lic,1).ne.0) pmax1=picblp
          if(nps(lic,2).ne.0) pmax2=picbls
          pmax=min(pmax1,pmax2)
          pmax=min(pmax,picbup)-ptol
          return
        else
c*** this for mantle-turning rays
          pmin1=0
          pmin2=0
          if(nps(loc-1,2).ne.0) pmin1=pcmbus
          if(nps(loc-1,1).ne.0) pmin2=pcmbup
          pmin=max(pmin1,pmin2)+ptol
c          if(iup.ne.0)then
cc  p or s at source
c            y=rdep/rnorm
c            call getmod(y,lsrc,rho,vpv,vsv,vph,vsh,eta)
c            if(iup.eq.1) pmax=rdep/min(vpv,vph)
c            if(iup.eq.2) pmax=rdep/min(vsv,vsh)
c          else 
c  P or S at source
c*** make lower-mantle turning only
            pmax1=1.e10
            pmax2=1.e10
            if(nps(lsrc,1).ne.0) pmax1=p660lp
            if(nps(lsrc,2).ne.0) pmax2=p660ls
            pmax=min(pmax1,pmax2)-ptol
c          endif
        end if
      end if
      return                      
C *** error in ray name
   20 ierr=1
      return
      END
 
      function ichdec(chr,ktype)
c Decodes characters in ray name, returning an integer.
      character*1 chr
      character*10 okchr                                              
      okchr(1:10)=' PSIJKpsic'
      DO 1 I=1,10
      ichdec=i
      if(chr.eq.okchr(i:i)) goto 2
    1 continue
      ichdec=11                 
      return
    2 ktype=1
      if(ichdec.eq.3.or.ichdec.eq.5.or.ichdec.eq.8) ktype=2
      return
      end
 
      subroutine addray (nps,ic1,ic2,i)
c Keeps sum of number of P ans S legs in each layer of the earth model.
      dimension nps(*)
      do 1 j=ic1,ic2
    1 nps(j)=nps(j)+i
      return
      end 

      subroutine zero(z,a1,b1,f,ftarg)
c  finds a root (a=z) of  f(a)-ftarg=0 with z
c  between a1,b1 by a combination of bisection and lin interp.
      implicit double precision(a-h,o-z)
	common/checkzero/nto
      data re/1.d-10/
      a=a1
      fa=f(a)-ftarg
      if(fa.eq.0.d0) then
         z=a
         return
      end if
      b=b1
      fb=f(b)-ftarg
      if(fb.eq.0.d0) then
         z=b
         return
      end if
c*** return if no zero (or not monotonic??)
      if(fa*fb.ge.0.d0) then
        z=-1.d0
        return
      end if
      c=a
      fc=fa
      s=c
      fs=fc

c===============================LAPO============================
	NLAPO=0
	NTO=0
c===============================LAPO============================

c ==============================================================
   10 h=0.5d0*(b+c)
	
c===============================LAPO============================
	NLAPO=NLAPO+1
	IF(NLAPO.GT.1000)THEN
	   NTO=1
	   RETURN
	ENDIF
c===============================LAPO============================

      t=dabs(h*re)
      if(dabs(h-b).le.t) then
        z=h
        return
      end if 
      if(dabs(fb).gt.dabs(fc)) then
        y=b
        fy=fb
        g=b
        fg=fb
        s=c
        fs=fc
      else
        y=s
        fy=fs
        g=c
        fg=fc
        s=b
        fs=fb
      end if
      if(fy.eq.fs) then
        b=h
      else
        e=(s*fy-y*fs)/(fy-fs)
        if(dabs(e-s).le.t) e=s+dsign(t,g-s)
        if((e-h)*(s-e).lt.0.d0) then
          b=h
        else
          b=e
        end if
      end if
      fb=f(b)-ftarg
      if(fg*fb.ge.0.d0) then
        c=s
        fc=fs
      else
        c=g
        fc=fg
      end if
      goto 10
c =============================================================
      end

      subroutine zero1(z,a1,b1,iq,ierr)
c  'z' is the zero of qtau (ie turning radius)  placed between 'a1'
c  and 'b1'. it is found by repeated linear interpolation between the
c  bounds of smaller and smaller intervals around the zero.
      implicit double precision(a-h,o-z)
      data re/1.d-14/
      a=a1
      call qtau(a,iq,fa)
      if(fa.eq.0.d0) then
         z=a
         return
      end if
      b=b1
      call qtau(b,iq,fb)
      if(fb.eq.0.d0) then
         z=b
         return
      end if
c******* no zero crossing or not monotonic
      if(fa*fb.ge.0.d0) then
         z=0.d0
         ierr=1
         return
      end if
      ierr=0
      c=a
      fc=fa
      s=c
      fs=fc
c =======================================================
c--this block to be repeated until a satisfactory estimate is reached
   10 h=0.5d0*(b+c)
      t=dabs(h*re)
      if(dabs(h-b).le.t) then
        z=h
        return
      end if
      if(dabs(fb).gt.dabs(fc)) then
        y=b
        fy=fb
        g=b
        fg=fb
        s=c
        fs=fc
      else
        y=s
        fy=fs
        g=c
        fg=fc
        s=b
        fs=fb
      end if
      if(fy.eq.fs) then
        b=h
      else 
        e=(s*fy-y*fs)/(fy-fs)
        if(dabs(e-s).le.t) e=s+dsign(t,g-s)
        if((e-h)*(s-e).lt.0.d0) then
          b=h
        else
          b=e
        end if
      end if
      call qtau(b,iq,fb)
      if(fg*fb.ge.0.d0) then
        c=s
        fc=fs
      else
        c=g
        fc=fg
      end if
      goto 10
c ================================================================
      end

      subroutine raypath(ray)
c Reads the ray name, breaks the ray down into its constituent parts, and calls
c subroutine addpath to compute the ray path. No checking of validity of ray
c name is performed as this has been done in raynam.
      implicit real*8(a-h,o-z)
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/path$/delt(2,800,2),nfin(2),kdep
      common/slow/picblp,picbls,picbup,pcmblp,pcmbup,pcmbus,
     +   p660lp,p660ls,lic,loc,llm,lsrc,klay(20)
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay

c---------------------------------LAPO----------------------------------
	common/timetot/timetot,kdown,kupn,iswi
	common/cmbradius/rcmbkm
	common/anikertot/dqdvphtot,dqdvpvtot,dqdvshtot,
     &dqdvsvtot,dqdetatot
c---------------------------------LAPO----------------------------------

      character*(*) ray      
      data if1/1/

c---------------------------------LAPO----------------------------------
	rcmbkm=xt(loc)
c---------------------------------LAPO----------------------------------

c*** write out circles if this is first call
      if(if1.eq.1) then
        if1=0
        ricb=xt(lic)/6371.d0
C        call circle(ricb)
        rmcb=xt(loc)/6371.d0
C        call circle(rmcb)
        rsur=1.
C        call circle(rsur)
      end if
c*** track ray legs
      th=0.d0

c---------------------------------LAPO----------------------------------
	timetot=0.d0
	dqdvphtot=0.d0
	dqdvpvtot=0.d0
	dqdvshtot=0.d0
	dqdvsvtot=0.d0
	dqdetatot=0.d0
	iswi=0
	kdown=0
	kupn=0
c---------------------------------LAPO----------------------------------

      j1=0
      do 10 i=1,70
      j=ichdec(ray(i:i),k)
      goto(10,11,11,12,12,13,14,14,100,100,10),j
C**** P or S ****
   11 goto(10,111,112,10,10,114,113,113,10,100),j1
      call addpath(th,k,kdep,klay(loc-1),1)
      goto 100
  111 call addpath(th,1,1,klay(loc-1),-1)
      call addpath(th,k,1,klay(loc-1),1)
      goto 100
  112 call addpath(th,2,1,klay(loc-1),-1)
  113 call addpath(th,k,1,klay(loc-1),1)
      goto 100
  114 call addpath(th,1,klay(loc-1)+1,klay(lic-1),-1)
      goto 100
C**** I or J ****
   12 call addpath(th,k,klay(lic-1)+1,klay(numlyr),1)
      call addpath(th,k,klay(lic-1)+1,klay(numlyr),-1)
      goto 100
C**** K ****
   13 goto(20,132,132,100,100,131,20,20,100,20),j1
  131 call addpath(th,1,klay(loc-1)+1,klay(lic-1),-1)
  132 call addpath(th,1,klay(loc-1)+1,klay(lic-1),1)
      goto 100
C**** p or s ****
   14 call addpath(th,k,1,kdep,-1)
  100 j1=j
      k1=k
   10 continue
C****  receiver leg ****
      call addpath(th,k1,1,klay(loc-1),-1)
c	write(10,*) '9999   9999'
   20 return                  
      end
 
c=======================================================================
      subroutine addpath(th,k,i1,i2,kup)
c  input:
c   k       = 1 for P leg, 2 for S leg
c   i1,i2 = indices in depth to go between
c   kup     = 1 for downgoing wave, -1 for upgoing
c   th is the current distance of raypath
      implicit real*8(a-h,o-z)
      common/path$/delt(2,800,2),nfin(2),kdep

c---------------------------------LAPO----------------------------------
c--I checked (10/98): the sum of the anisotropic integrated vph
c--and vpv kernels seems to be equal to the total travel time array.
c--that is, at any th, timetot=-(dqdvphtot+dqdvpvtot).

	COMMON/TIMEINCR/TIMEARRAY(800,2)

	COMMON/ANKERINC/DQDVPH(800,2),DQDVPV(800,2),DQDVSH(800,2),
     &DQDVSV(800,2),DQDETA(800,2)
	common/anikertot/dqdvphtot,dqdvpvtot,dqdvshtot,
     &dqdvsvtot,dqdetatot

c--timetot is the current cumulative travel time
	common/timetot/timetot,kdown,kupn,iswi
	parameter(ndimmax2=10000)
	common/lapo/raddown(ndimmax2),radup(ndimmax2),kdowntot,kuptot,cultime(ndimmax2),
     &dedown(ndimmax2),deup(ndimmax2),culdel(ndimmax2),culray(ndimmax2),iswtc

	common/lapoanis/CULdqdvph(ndimmax2),CULdqdvpv(ndimmax2),CULdqdvsh(ndimmax2),
     &CULdqdvsv(ndimmax2),CULdqdeta(ndimmax2)

c--should work for core phases, not for phases with multiple legs.
	data pi/3.14159265358979d0/
c---------------------------------LAPO----------------------------------

      i2s = min(i2,nfin(k))
      if (kup.eq.1) then

c downgoing leg
        do l=i1,i2s
          rray=delt(1,l,k)/6371.
          th=th+delt(2,l,k)

c---------------------------------LAPO----------------------------------
	  timetot=timetot+timearray(l,k)
	  dqdvphtot=dqdvphtot+DQDVPH(l,k)
	  dqdvpvtot=dqdvpvtot+DQDVPV(l,k)
	  dqdvshtot=dqdvshtot+DQDVsH(l,k)
	  dqdvsvtot=dqdvsvtot+DQDVsV(l,k)
	  dqdetatot=dqdetatot+DQDeta(l,k)

          if(rray.ne.0.d0) then
	     KDOWN=KDOWN+1
	     RADDOWN(KDOWN)=DELT(1,L,K)
	     CULTIME(KDOWN)=TIMETOT
	     CULdqdvph(KDOWN)=dqdvphtot
	     CULdqdvpv(KDOWN)=dqdvpvtot
c	write(46,*)th*180./pi,timetot
c	write(47,*)th*180./pi,-(dqdvphtot+dqdvpvtot)
c	write(48,*)th*180./pi,-(dqdvpvtot),-dqdvphtot,timetot
C	write(10,*)rray*sin(th),rray*cos(th)
	     CULdqdvsh(KDOWN)=dqdvshtot
	     CULdqdvsv(KDOWN)=dqdvsvtot
	     CULdqdeta(KDOWN)=dqdetatot
	     CULray(KDOWN)=DELT(1,L,K)
	     CULdel(KDOWN)=th
	     dedown(KDOWN)=th
	  endif
c---------------------------------LAPO----------------------------------

	enddo

c---------------------------------LAPO----------------------------------
c--radup and cultime will not be consistent but that shouldnt matter. 
c--notice that the following happens only if the ray has a 
c--downgoing leg. all this will not work for ray w/ multiple legs.
	iswi=1
	radup(1)=RADDOWN(KDOWN)
	deup(1)=dedown(kdown)


c--(the "last" value of kdown counts)
	kdowntot=kdown
c---------------------------------LAPO----------------------------------

      else
c upgoing leg
        i1s=max(i1,2)
        do l=i2s,i1s,-1
          rray=delt(1,l-1,k)/6371.
          th=th+delt(2,l,k)

c---------------------------------LAPO----------------------------------
	  timetot=timetot+timearray(l,k)
	  dqdvphtot=dqdvphtot+DQDVPH(l,k)
	  dqdvpvtot=dqdvpvtot+DQDVPV(l,k)
	  dqdvshtot=dqdvshtot+DQDVsH(l,k)
	  dqdvsvtot=dqdvsvtot+DQDVsV(l,k)
	  dqdetatot=dqdetatot+DQDeta(l,k)

          if(rray.ne.0.d0)then
	     Kupn=Kupn+1
	     RADup(Kupn+iswi)=DELT(1,L-1,K)
	     CULTIME(KDOWN+kupn)=TIMETOT
	     CULdqdvph(KDOWN+kupn)=dqdvphtot
	     CULdqdvpv(KDOWN+kupn)=dqdvpvtot
c	write(46,*)th*180./pi,timetot
c	write(47,*)th*180./pi,-(dqdvphtot+dqdvpvtot)
c	write(48,*)th*180./pi,-(dqdvpvtot),-dqdvphtot,timetot
C	write(10,*)rray*sin(th),rray*cos(th)
	     CULdqdvsh(KDOWN+kupn)=dqdvshtot
	     CULdqdvsv(KDOWN+kupn)=dqdvsvtot
	     CULdqdeta(KDOWN+kupn)=dqdetatot
	     CULray(KDOWN+kupn)=DELT(1,L-1,K)
	     CULdel(KDOWN+kupn)=th
	     deup(kupn+iswi)=th
	  endif
c---------------------------------LAPO----------------------------------

        enddo

c---------------------------------LAPO----------------------------------
	kuptot=kupn
	iswtc=iswi
c---------------------------------LAPO----------------------------------

      endif
      return
      end
c=======================================================================

      subroutine circle(radi)
      implicit real*8(a-h,o-z)
      data pi/3.14159265358979d0/
      rad=180./pi
      do i=1,361
        ang=i/rad
C        write(10,*) radi*sin(ang),radi*cos(ang)
      enddo
      return
      end
