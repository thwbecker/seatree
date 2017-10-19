	include 'common_para.h'
c	parameter(eq_incr=3)
c	parameter(ifa=1)
c	parameter(n=4592)
	parameter(n=15000) ! maximum number of free parameters

c	parameter(eq_hi=eq_incr/ifa)
c	parameter(m=100000,nonz=100000000)
	parameter(m=500000,nonz=200000000)
c	parameter(nlatzones=180./eq_incr)
c	parameter(nlatzohi=180./eq_hi)
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
	dimension nsqrsh(nlatzhmax),nsqtoth(nlatzhmax+1)
	parameter(radian=pi/180.0)
	dimension v(n),w(n),x(n),se(n),aty(n)
	dimension t(m),ttt(m)
	dimension indx(nonz),values(nonz),mpoin(0:m)
	character answ*1
c	dimension nsqrs(nlatzones),nsqtot(nlatzones+1)
c	dimension nsqrsh(nlatzohi),nsqtoth(nlatzohi+1)
	dimension iresolu(500000),ire2(500000)
	dimension inew(500000),iold(500000),inewh(500000),ioldh(500000)
	character type*1,wave_type*1,perio4*4,chn*4
	character chincr*2,chfact*2,perio*3
	character filen*200,txtfile*200,fileout*200
	character dampchar1*5,dampchar2*5
       
	print*,"original data file?(use quotes if specifying whole path)"
	read*,txtfile
	print *,'using ',txtfile

	print*,"root of name of matrix files?"
	read(*,'(a80)')filen
	print *,'using ',filen

	print*,"solution file?"
	read(*,'(a80)')fileout
	print *,'using ',fileout

	print*,"size of coarse pixels?"
	read*,eq_incr
	print *,'using ',eq_incr
	
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

	print*,'number of latitudinal zones: ',nlatzones
	if(mod(nlatzones,2).ne.0)then
	   print*,'it should not be odd'
	   stop
	endif

c--------------------------determine the arrays that the define the grid:
	call parametra2(eq_incr,nsqrs,nsqtot,nsqrsh,nsqtoth,
     &	nlatzones,eq_hi,nlatzohi,ifa,iswit,refgrid)
c--determine number of coarse blocks and number of fine blocks:
	nuto=nsqtot(nlatzones+1)
	nuhi=nsqtoth(nlatzohi+1)
	call numbers(nco,nfi,westbo,eastbo,southbo,rthnobo,
     &	nlatzones,nlatzohi,nsqrs,nsqtot,nsqrsh,nsqtoth,eq_incr,
     &	eq_hi,nuto,nuhi,ifa,iresolu,kireso,ire2,kire2)
	print*,'in the actual grid there are'
	print*,nco,' coarse blocks'
	print*,nfi,' fine blocks'
	print*,nco+nfi,' total blocks'
c	if(nco+nfi.ne.n)then
c	   print*,'but n=',n
c	   pause
c	endif
	ieqi=eq_incr
	write(chn,"(i2.2)")ieqi
c	open(99,file="n_"//chn)
	open(99,file="nome")
	write(99,*)nco+nfi
	close(99)
c--determine inew,inewh,iold,ioldh
	call correspc(iresolu,kireso,nuto,0,inew,iold,ico)
	print*,'compare:',nco,inew(nuto)
	iof=nco
	call corresph(ire2,kire2,nuhi,iof,inewh,ioldh,ico)
c--remove the average
	print*,'remove the average from the solution?'
	read*,answ
c--choose norm damping parameters:
	print*,'norm minimization parameter?'
	read*,damp0
c--correct damp for consistency with other parametrizations
c--(boschi and dziewonski, jgr 1999):
	omega=eq_incr*radian
	damp=damp0*omega
c	damp=damp0	! not anymore

c--in the high-resolution region we should use a different correction
	omega1=eq_hi*radian
c	damp1=damp0*omega1
	damp1=damp0*(omega1/omega)
c--choose roughness damping parameter:
	print*,'roughness minimization parameter?'
	read*,gradamp
	print*,'roughness damping parameter within the
     & high resolution region?'
	read*,gradamp_fin

c-------------find reference velocity (look at original data file from erik)
c	write(perio4,"(i4.4)")input_per
c	txtfile=
c     &'/home/simona/tomografia/etl_summary_phase_data/wei_sum.02.'
c     &//type//perio4//'_1.txt'

	open(2,file=txtfile,status='old')
	read(2,*)wave_type,period,velo0
c--qui sopra usa format 1003 se cosi non funziona!
	print*,wave_type,period,velo0
	int_period=period
	write(perio4,"(i4.4)")int_period
	close(2)
1003	format(8x,a1,8x,f7.4,10x,f8.5)

	print*,'norm damping parameters=',damp,' and (hi-res)',damp1
	print*,'roughness damping parameters=',gradamp,gradamp_fin

c-----------------------------read in the matrix

c	ieq_incr=eq_incr
c	write(chincr,'(i2.2)')ieq_incr
c	write(chfact,'(i2.2)')ifa
c	filen=type//'_'//perio//'_'//chincr//'_'//chfact

	do k=1,80
	   if(filen(k:k).eq." ")goto 74
	enddo

74	kend=k-1
	open(1,file=filen(1:kend)//'.xxx',status='old',
     &   access='direct',recl=4,
     &  form='unformatted')
	open(4,file=filen(1:kend)//'.ind',status='old',
     &   access='direct',recl=4,
     &  form='unformatted')
	open(77,file=filen(1:kend)//'.rhs',status='old')
	open(3,file=filen(1:kend)//'.pnt',status='old')

c--read matrix
	mpoin(0)=0
	do icol=1,m
	   read(3,*,end=9)mpoin(icol)
	   read(77,*,end=7)t(icol)
	   do jj=mpoin(icol-1)+1,mpoin(icol)
	      read(4,rec=jj)indx(jj)
	      read(1,rec=jj)values(jj)
	   enddo
	enddo
7	print*,'rhs was shorter or m too small !!!'
9	print*,'number of rows=',icol-1
	kp=icol-1
	kp0=kp
	nz=jj-1
	close(1)
	close(4)
	close(3)
	close(77)

c--save r.h.s.; we are going to need it to compute variance reduction.
	print*,'save r.h.s.'
	do it=1,kp0
	   ttt(it)=t(it)
	enddo

c--impose norm minimization of coarser grid (stronger for consistency)
	nelp=kp
	nnn=nz
	dampco=damp-damp1
	print*,'coarse grid norm damping parameter=',dampco
c--loop over all coarse blocks:
	do icoblo=1,nco
	   nnn=nnn+1
	   values(nnn)=dampco
	   indx(nnn)=icoblo
	   nnn=nnn+1
	   nelp=nelp+1
	   mpoin(nelp)=nnn
	   t(nelp)=0.
37	continue
	enddo
	nz=nnn
	kp=nelp
	kp1=kp

c--impose a "total gradient" (roughness) minimization constraint
	if(gradamp.eq.0.)goto13
c--damp low-resolution part of the model
	call rough_d_outside(values,indx,mpoin,n,gradamp,kp,nz,nonz,
     &	m,t,nsqrs,nsqtot,nlatzones,eq_incr,n,iresolu,kireso,inew,
     &	eq_hi,ifa,nsqrsh,nsqtoth,nlatzohi,inewh,nuhi)
c--damp high-resolution part of the model
	call rough_d_inside(values,indx,mpoin,n,gradamp_fin,kp,nz,nonz,
     &	m,t,nsqrsh,nsqtoth,nlatzohi,eq_hi,nuhi,
     &	nsqrs,nsqtot,nlatzones,inew,inewh,ifa,eq_incr)

c--internal check for consistency
	do i=kp1+1,kp
c	   print*,'row:',i
	   j1=mpoin(i-1)+1
	   j2=mpoin(i)
	   sumrow=0.
	   do j=j1,j2
	      sumrow=sumrow+values(j)
	   enddo
	   if(sumrow.gt.0.0001)then
	      print*,'error in damping matrix, row',i,' sum=',sumrow
	      stop
	   endif
c	   write(*,*)(indx(j),j=j1,j2)
c	   write(*,*)(values(j),j=j1,j2)
	enddo
c	pause

13	atol=0.
	btol=0.
	conlim=0.
	itnlim=10*(nco+nfi)
	nout=1

	npx=nco+nfi

	print*,'general norm damping parameter=',damp1
	print*,'call lsqr'
	call lsqr( m, npx, damp1,
     $             1, npx, iw, aty,
     $             t, v, w, x, se,
     $             atol, btol, conlim, itnlim, nout,
     $       istop, itn, anorm, acond, rnorm, arnorm, xnorm,
     &       indx,values,mpoin,nonz)
cTEST
c	print*,"istop=",istop
	if(istop.ne.5)stop "does not converge"

	czero=velo0/(6371.*radian)
c--remove the average
	ave=0.
c	if((answ.eq.'y').or.(answ.eq.'Y'))then
c	   print*,'remove the average'
	   do jj=1,npx
	      ave=ave+x(jj)
	   enddo
	   ave=ave/npx
	   print*,'average=',(-czero*ave*100.),czero,ave
    	   write(dampchar1,'(i5.5)')int(damp0)
	   write(dampchar2,'(i5.5)')int(gradamp)
c	   open(72,file="averages."//dampchar1//"_"//dampchar2//".txt",access="append")
	   open(72,file="averages.txt",access="append")
c	   write(72,*)"------"
c	   write(72,*)period,(-czero*ave*100.)," ",filen," ",damp0,gradamp
	write(72,*)damp0,gradamp,-czero*ave*100.,itn,eq_incr
	   close(72)
c	else
c	   print*,'keep the average'
c	endif
	if(answ.ne."y".and.answ.ne."Y")ave=0.

c--compute variance reduction
	print*,'compute variance reduction'
	rnumer=0.
	denom=0.
	do l=1,kp0
	   denom=denom+(ttt(l)*ttt(l))
	   tot=0.
	   do ll=mpoin(l-1)+1,mpoin(l)
	      tot=tot+(values(ll)*x(indx(ll)) )
	   enddo
	   rnumer=rnumer+(ttt(l)-tot)**2
	enddo
	open(72,file="loglike.txt",access="append")
c--log-likelihood function
	write(72,*)damp0,gradamp,-0.5*float(kp0)*log(rnumer),itn,npx
	close(72)
	varred=1.-(rnumer/denom)
	write(*,*)'variance reduction=',varred,' norm=',xnorm
c	open(72,file="fit."//dampchar1//"_"//dampchar2//".txt",access="append")
c	open(72,file="fit_"//wave_type//"_"//perio4//".txt",access="append")
	open(72,file="fit.txt",access="append")
c	write(72,*)filen
	write(72,*)damp0,gradamp,varred,itn,eq_incr
c	write(72,*)"---"
	close(72)

c	write(dampchar1,'(i5.5)')int(damp0)
c	write(dampchar2,'(i5.5)')int(gradamp)
	open(1,file=fileout)
	xnorm_my=0.0d0
	do jj=1,npx
c	result=1./(1./czero+(x(jj)-ave))-czero
c	result=result/czero
c	   result=result*100.
c	   write(1,*)jj,result

	   write(1,*)jj,(x(jj)-ave)*100.
	   xnorm_my=xnorm_my+x(jj)**2

	enddo
	xnorm_my=sqrt(xnorm_my)
	write(1,*)'iter.s:',itn,' v.red.:',varred,' damp:',
     &damp0,gradamp,gradamp_fin
c	print *,'norm ',xnorm,xnorm_my
	close(1)

c=====================================================================
	stop
	end


c================================================================
c================================================================
* from rfischer@seismology.harvard.edu wed jul 23 15:30:22 1997
* retrieved by bob fischer from: http://www.netlib.org/linalg/lsqr
* from arpa!sol-michael.stanford.edu!mike 5 may 89 23:53:00 pdt
      subroutine lsqr  ( m, n, damp,
     $                   leniw, lenrw, iw, rw,
     $                   u, v, w, x, se,
     $                   atol, btol, conlim, itnlim, nout,
     $   istop, itn, anorm, acond, rnorm, arnorm, xnorm,
     &       indx,values,mpoin,nonz)

c--i had to add next line to pass the matrix to aprod without using a common
	dimension indx(nonz),values(nonz),mpoin(0:m)

c      external           aprod
      integer            m, n, leniw, lenrw, itnlim, nout, istop, itn
      integer            iw(leniw)
c      double precision   rw(lenrw), u(m), v(n), w(n), x(n), se(n),
      real*4  rw(lenrw), u(m), v(n), w(n), x(n), se(n),
     $                   atol, btol, conlim, damp,
     $                   anorm, acond, rnorm, arnorm, xnorm

*-----------------------------------------------------------------------
*     intrinsics and local variables

      intrinsic          abs, mod, sqrt
      integer            i, nconv, nstop
c      double precision   dnrm2
	real*4 dnrm2
c      double precision   alfa, bbnorm, beta, bnorm,
      real*4   alfa, bbnorm, beta, bnorm,
     $                   cs, cs1, cs2, ctol, dampsq, ddnorm, delta,
     $                   gamma, gambar, phi, phibar, psi,
     $                   res1, res2, rho, rhobar, rhbar1, rhbar2,
     $                   rhs, rtol, sn, sn1, sn2,
     $                   t, tau, test1, test2, test3,
     $                   theta, t1, t2, t3, xxnorm, z, zbar

c      double precision   zero,           one
      parameter        ( zero = 0.,  one = 1. )

      character*16       enter, exit
      character*60       msg(0:7)

      data               enter /' enter lsqr.    '/,
     $                   exit  /' exit  lsqr.    '/

      data               msg
     $ / 'the exact solution is  x = 0',
     $   'ax - b is small enough, given atol, btol',
     $   'the least-squares solution is good enough, given atol',
     $   'the estimate of cond(abar) has exceeded conlim',
     $   'ax - b is small enough for this machine',
     $   'the least-squares solution is good enough for this machine',
     $   'cond(abar) seems to be too large for this machine',
     $   'the iteration limit has been reached' /
*-----------------------------------------------------------------------


*     initialize.

      if (nout .gt. 0)
     $   write(nout, 1000) enter, m, n, damp, atol, conlim, btol, itnlim
      itn    =   0
      istop  =   0
      nstop  =   0
      ctol   =   zero
      if (conlim .gt. zero) ctol = one / conlim
      anorm  =   zero
      acond  =   zero
      bbnorm =   zero
      dampsq =   damp**2
      ddnorm =   zero
      res2   =   zero
      xnorm  =   zero
      xxnorm =   zero
      cs2    = - one
      sn2    =   zero
      z      =   zero

      do 10  i = 1, n
         v(i)  =  zero
         x(i)  =  zero
        se(i)  =  zero
   10 continue

*     set up the first vectors u and v for the bidiagonalization.
*     these satisfy  beta*u = b,  alfa*v = a(transpose)*u.

      alfa   =   zero
      beta   =   dnrm2 ( m, u, 1 )

      if (beta .gt. zero) then
         call dscal ( m, (one / beta), u, 1 )
         call aprod ( 2, m, n, v, u, leniw, lenrw, iw, rw ,
     &       indx,values,mpoin,nonz)
         alfa   =   dnrm2 ( n, v, 1 )
      end if

      if (alfa .gt. zero) then
         call dscal ( n, (one / alfa), v, 1 )
         call dcopy ( n, v, 1, w, 1 )
      end if

      arnorm =   alfa * beta
      if (arnorm .eq. zero) go to 800

      rhobar =   alfa
      phibar =   beta
      bnorm  =   beta
      rnorm  =   beta

      if (nout   .gt.  0  ) then
         if (dampsq .eq. zero) then
             write(nout, 1200)
         else
             write(nout, 1300)
         end if
         test1  = one
         test2  = alfa / beta
         write(nout, 1500) itn, x(1), rnorm, test1, test2
         write(nout, 1600)
      end if

*     ------------------------------------------------------------------
*     main iteration loop.
*     ------------------------------------------------------------------
  100 itn    = itn + 1
	print*,'iteration:',itn

*     perform the next step of the bidiagonalization to obtain the
*     next  beta, u, alfa, v.  these satisfy the relations
*                beta*u  =  a*v  -  alfa*u,
*                alfa*v  =  a(transpose)*u  -  beta*v.

      call dscal ( m, (- alfa), u, 1 )
      call aprod ( 1, m, n, v, u, leniw, lenrw, iw, rw,
     &       indx,values,mpoin,nonz)

      beta   =   dnrm2 ( m, u, 1 )
      bbnorm =   bbnorm  +  alfa**2  +  beta**2  +  dampsq

      if (beta .gt. zero) then
         call dscal ( m, (one / beta), u, 1 )
         call dscal ( n, (- beta), v, 1 )
         call aprod ( 2, m, n, v, u, leniw, lenrw, iw, rw ,
     &       indx,values,mpoin,nonz)
         alfa   =   dnrm2 ( n, v, 1 )
         if (alfa .gt. zero) then
            call dscal ( n, (one / alfa), v, 1 )
         end if
      end if

*     use a plane rotation to eliminate the damping parameter.
*     this alters the diagonal (rhobar) of the lower-bidiagonal matrix.

      rhbar2 = rhobar**2  +  dampsq
      rhbar1 = sqrt( rhbar2 )
      cs1    = rhobar / rhbar1
      sn1    = damp   / rhbar1
      psi    = sn1 * phibar
      phibar = cs1 * phibar

*     use a plane rotation to eliminate the subdiagonal element (beta)
*     of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.

      rho    =   sqrt( rhbar2  +  beta**2 )
      cs     =   rhbar1 / rho
      sn     =   beta   / rho
      theta  =   sn * alfa
      rhobar = - cs * alfa
      phi    =   cs * phibar
      phibar =   sn * phibar
      tau    =   sn * phi

*     update  x, w  and the standard error estimates.

      t1     =   phi   / rho
      t2     = - theta / rho
      t3     =   one   / rho

      do 200  i =  1, n
         t      =  w(i)
         x(i)   =  t1*t  +  x(i)
         w(i)   =  t2*t  +  v(i)
         t      = (t3*t)**2
         se(i)  =  t     +  se(i)
         ddnorm =  t     +  ddnorm
  200 continue

*     use a plane rotation on the right to eliminate the
*     super-diagonal element (theta) of the upper-bidiagonal matrix.
*     then use the result to estimate  norm(x).

      delta  =   sn2 * rho
      gambar = - cs2 * rho
      rhs    =   phi    - delta * z
      zbar   =   rhs    / gambar
      xnorm  =   sqrt( xxnorm    + zbar **2 )
      gamma  =   sqrt( gambar**2 + theta**2 )
      cs2    =   gambar / gamma
      sn2    =   theta  / gamma
      z      =   rhs    / gamma
      xxnorm =   xxnorm + z**2

*     test for convergence.
*     first, estimate the norm and condition of the matrix  abar,
*     and the norms of  rbar  and  abar(transpose)*rbar.

      anorm  =   sqrt( bbnorm )
      acond  =   anorm * sqrt( ddnorm )
      res1   =   phibar**2
      res2   =   res2  +  psi**2
      rnorm  =   sqrt( res1 + res2 )
      arnorm =   alfa  * abs( tau )

*     now use these norms to estimate certain other quantities,
*     some of which will be small near a solution.

      test1  =   rnorm /  bnorm
      test2  =   zero
      if (rnorm .gt. zero) test2 = arnorm / (anorm * rnorm)
      test3  =   one   /  acond
      t1     =   test1 / (one  +  anorm * xnorm / bnorm)
      rtol   =   btol  +  atol *  anorm * xnorm / bnorm

*     the following tests guard against extremely small values of
*     atol, btol  or  ctol.  (the user may have set any or all of
*     the parameters  atol, btol, conlim  to zero.)
*     the effect is equivalent to the normal tests using
*     atol = relpr,  btol = relpr,  conlim = 1/relpr.

      t3     =   one + test3
      t2     =   one + test2
      t1     =   one + t1
      if (itn .ge. itnlim) istop = 7
      if (t3  .le. one   ) istop = 6
      if (t2  .le. one   ) istop = 5
      if (t1  .le. one   ) istop = 4

*     allow for tolerances set by the user.

      if (test3 .le. ctol) istop = 3
      if (test2 .le. atol) istop = 2
      if (test1 .le. rtol) istop = 1
*     ==================================================================

*     see if it is time to print something.

      if (nout  .le.  0       ) go to 600
      if (n     .le. 40       ) go to 400
      if (itn   .le. 10       ) go to 400
      if (itn   .ge. itnlim-10) go to 400
      if (mod(itn,10) .eq. 0  ) go to 400
      if (test3 .le.  2.0*ctol) go to 400
      if (test2 .le. 10.0*atol) go to 400
      if (test1 .le. 10.0*rtol) go to 400
      if (istop .ne.  0       ) go to 400
      go to 600

*     print a line for this iteration.

  400 write(nout, 1500) itn, x(1), rnorm, test1, test2, anorm, acond
      if (mod(itn,10) .eq. 0) write(nout, 1600)
*     ==================================================================

*     stop if appropriate.
*     the convergence criteria are required to be met on  nconv
*     consecutive iterations, where  nconv  is set below.
*     suggested value:  nconv = 1, 2  or  3.

  600 if (istop .eq. 0) nstop = 0
      if (istop .eq. 0) go to 100
      nconv  =   1
      nstop  =   nstop + 1
      if (nstop .lt. nconv  .and.  itn .lt. itnlim) istop = 0
      if (istop .eq. 0) go to 100
*     ------------------------------------------------------------------
*     end of iteration loop.
*     ------------------------------------------------------------------


*     finish off the standard error estimates.

      t    =   one
      if (m      .gt.   n )  t = m - n
      if (dampsq .gt. zero)  t = m
      t    =   rnorm / sqrt( t )

      do 700  i = 1, n
         se(i)  = t * sqrt( se(i) )
  700 continue

*     print the stopping condition.

  800 if (nout .gt. 0) then
         write(nout, 2000) exit, istop, itn,
     $                     exit, anorm, acond,
     $                     exit, rnorm, arnorm,
     $                     exit, bnorm, xnorm
         write(nout, 3000) exit, msg(istop)
      end if

  900 return

*     ------------------------------------------------------------------
 1000 format(// 1p, a, '  least-squares solution of  a*x = b'
     $    / ' the matrix  a  has', i7, ' rows   and', i7, ' columns'
     $    / ' the damping parameter is         damp   =', e10.2
     $    / ' atol   =', e10.2, 15x,        'conlim =', e10.2
     $    / ' btol   =', e10.2, 15x,        'itnlim =', i10)
 1200 format(// '   itn       x(1)           function',
     $   '     compatible   ls        norm a    cond a' /)
 1300 format(// '   itn       x(1)           function',
     $   '     compatible   ls     norm abar cond abar' /)
 1500 format(1p, i6, 2e17.9, 4e10.2)
 1600 format(1x)
 2000 format(/ 1p, a, 6x, 'istop =', i3,   16x, 'itn    =', i9
     $       /     a, 6x, 'anorm =', e13.5, 6x, 'acond  =', e13.5
     $       /     a, 6x, 'rnorm =', e13.5, 6x, 'arnorm =', e13.5,
     $       /     a, 6x, 'bnorm =', e13.5, 6x, 'xnorm  =', e13.5)
 3000 format( a, 6x, a )
*     ------------------------------------------------------------------
*     end of lsqr
      end


c======================================================================
c======================================================================

      subroutine  dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
c      double precision dx(*),dy(*)
      real dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end

c================================================================
c================================================================

      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
c      double precision da,dx(*)
      real*4 da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end

c===============================================================
c===============================================================

c      double precision function dnrm2 ( n, x, incx )
      real function dnrm2 ( n, x, incx )
*     .. scalar arguments ..
      integer                           incx, n
*     .. array arguments ..
c      double precision                  x( * )
      real*4               x( * )
*     ..
*
*  dnrm2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     dnrm2 := sqrt( x'*x )
*
*
*
*  -- this version written on 25-october-1982.
*     modified on 14-october-1993 to inline the call to dlassq.
*     sven hammarling, nag ltd.
*
*
*     .. parameters ..
c      double precision      one         , zero
      parameter           ( one = 1., zero = 0. )
*     .. local scalars ..
      integer               ix
c      double precision      absxi, norm, scale, ssq
      real*4     absxi, norm, scale, ssq
*     .. intrinsic functions ..
      intrinsic             abs, sqrt
*     ..
*     .. executable statements ..
      if( n.lt.1 .or. incx.lt.1 )then
         norm  = zero
      else if( n.eq.1 )then
         norm  = abs( x( 1 ) )
      else
         scale = zero
         ssq   = one
*        the following loop is equivalent to this call to the lapack
*        auxiliary routine:
*        call dlassq( n, x, incx, scale, ssq )
*
         do 10, ix = 1, 1 + ( n - 1 )*incx, incx
            if( x( ix ).ne.zero )then
               absxi = abs( x( ix ) )
               if( scale.lt.absxi )then
                  ssq   = one   + ssq*( scale/absxi )**2
                  scale = absxi
               else
                  ssq   = ssq   +     ( absxi/scale )**2
               end if
            end if
   10    continue
         norm  = scale * sqrt( ssq )
      end if
*
      dnrm2 = norm
      return
*
*     end of dnrm2.
*
      end

c================================================================
c================================================================

	subroutine aprod(mode,m,n,x,y,lenin,lenva,iw,aty,
     &       indx,values,mpoin,nonz)
c	implicit real*8(a-h,o-z)
	integer iw(lenin)
	real*4 x(n),y(m),aty(n)
	dimension indx(nonz),values(nonz),mpoin(0:m)

	if(mode.eq.1)then
c--compute y=y+a*x
	do k=1,m
	   pro=0.
	   do j=mpoin(k-1)+1,mpoin(k)
	      pro=pro+values(j)*x(indx(j))
	   enddo
	   y(k)=y(k)+pro
	enddo
	return

	elseif(mode.eq.2)then
	do i=1,n
	   aty(i)=0.
	enddo
c--compute x=x+(a^t)*y
	do k=1,m
	   do j=mpoin(k-1)+1,mpoin(k)
	      aty(indx(j))=aty(indx(j))+values(j)*y(k)
	   enddo
	enddo
	do i=1,n
	   x(i)=x(i)+aty(i)
	enddo
	return

	else
	print*,'error: mode=',mode
	stop
	endif
	end

c*************
	subroutine rough_d_outside(g,indx,mpoin,nblo,weight,kpoin,nnn,
     &	nonz,m,rhs,nsqrs,nsqtot,nlatzones,eq_incr,n,iresolu,kireso,inew,
     &	eq_hi,ifa,nsqrsh,nsqtoth,nlatzohi,inewh,numhi)
c--defines the matrix corresponding to gradient damping,
c--accounting for the surface gradient at each boundary between blocks.
	dimension indx(nonz),g(nonz),mpoin(0:m),rhs(m)
	dimension iresolu(10000),inew(50000)
c	dimension nsqrs(nlatzones),nsqtot(nlatzones+1)
c	dimension nsqrsh(nlatzohi),nsqtoth(nlatzohi+1)
	parameter(nlatzomax=180,nlatzhmax=720)
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
	dimension nsqrsh(nlatzhmax),nsqtoth(nlatzhmax+1)
	dimension inewh(50000)

	kpoin0=kpoin
	print*,'gradamp=',weight

	nnorth=n/2
c	print*,'nnorth=',nnorth
	indexch=(nlatzones/2)+1
	nsouth=nnorth+nsqrs(indexch)+1
c	print*,'nsouth=',nsouth

	do i1=1,nlatzones
c--rloin_2=half the longitudinal size of the block
	   rloin_2=180./nsqrs(i1)
	   rloinhi_2=rloin_2/ifa
	   ifirst=nsqtot(i1)+1
	   ilast=nsqtot(i1+1)
	   do i2=ifirst,ilast
	      do k=1,kireso
	         if(i2.eq.iresolu(k))then
	            goto47
	         endif
	      enddo
c--find iright
	      iright=i2+1
	      if(i2.eq.ilast)iright=ifirst
	      do k=1,kireso
	         if(iright.eq.iresolu(k))then
c--block iright is in the high-resolution area; 
c--find the corresponding easternmost finer blocks:
c	print*,i2,iright,inew(i2),inew(iright)
	            call coordsuper(i2,rila,rilo,nsqrs,nlatzones,eq_incr)
c	            ilorihi=rilo*100.+(rloin_2+rloinhi_2)*100.
c		    rloin=rilo*100.+(rloin_2+rloinhi_2)
		    rloin=rilo+(rloin_2+rloinhi_2) ! fixed 2/2006
c	            ilariastart=(rila+(eq_incr-eq_hi)/2.)*100.
		    rlainstart=rila+(eq_incr-eq_hi)/2.
	            do kuks=1,ifa
c	               ilaria=ilariastart-((kuks-1)*eq_hi)*100.
		       rlain=rlainstart-((kuks-1)*eq_hi)
c	               irihi=superisqre(ilaria,ilorihi,
c     &	nsqrsh,nsqtoth,nlatzohi,numhi,eq_hi)
	               irihi=isqre(rlain,rloin,
     &                   nsqrsh,nsqtoth,
     &                   nlatzohi,numhi,eq_hi)
c	print*,irihi,inewh(irihi)
	               indx(nnn+1)=inewh(irihi)
	               g(nnn+1)=-weight/ifa
	               nnn=nnn+1
	            enddo
c	pause
	            goto37
	         endif
	      enddo

c--define the row of the damping matrix corresponding to longitude variation
	      indx(nnn+1)=inew(iright)
	      g(nnn+1)=-weight
	      nnn=nnn+1
37	      indx(nnn+1)=inew(i2)
	      g(nnn+1)=weight
	      nnn=nnn+1
	      kpoin=kpoin+1
	      mpoin(kpoin)=nnn
	      rhs(kpoin)=0.

c==============================================
c==============================================
c--find iup or idw
	      if((1.le.i2).and.(i2.le.nnorth))then
c--if i2 is in the n hemisphere
	         call coordsuper(i2,rila,rilo,nsqrs,nlatzones,eq_incr)
	         dwla=rila-eq_incr
	         rilol=rilo-rloin_2
	         rilor=rilo+rloin_2
	         idwla=dwla*100.
	         idwl=isqre(dwla,rilol,nsqrs,nsqtot,nlatzones,n,eq_incr)
	         idwr=isqre(dwla,rilor,nsqrs,nsqtot,nlatzones,n,eq_incr)
c--define the row corresponding to
c--variation wrt latitude (n hemisphere)
	         g(nnn+1)=weight
	         indx(nnn+1)=inew(i2)
	         nnn=nnn+1
	         peso=1./(idwr-idwl+1)
	         do idw=idwl,idwr
	            do k=1,kireso
	               if(idw.eq.iresolu(k))then
	                  call rang(idw,xlamin,xlamax,xlomin,xlomax,
     &	nsqrs,nsqtot,nlatzones,n,eq_incr)
c--block idw is in the high-resolution area; 
c--find the corresponding northernmost finer blocks:
	                  fin_incr=(360./nsqrs(i1+1))/ifa
c	                  iladwhi=rila*100.-((eq_incr+eq_hi)/2.)*100.
			  rladwhi=rila-(eq_incr+eq_hi)/2.
c	                  ilostart=(xlomin+(fin_incr/2.))*100.
			  rlostart=xlomin+(fin_incr/2.)
	                  do kuks=1,ifa
c	                     ilodwhi=ilostart+((kuks-1)*(fin_incr))*100.
	                     rlodwhi=rlostart+(kuks-1)*(fin_incr)
c	                     idwhi=superisqre(iladwhi,ilodwhi,
c     &	nsqrsh,nsqtoth,nlatzohi,numhi,eq_hi)
	                     idwhi=isqre(rladwhi,rlodwhi,
     &	nsqrsh,nsqtoth,nlatzohi,numhi,eq_hi)
	                     indx(nnn+1)=inewh(idwhi)
	                     g(nnn+1)=-weight*peso/ifa
	                     nnn=nnn+1
	                  enddo
	                  goto57
	               endif
	            enddo
	            g(nnn+1)=-weight*peso
	            indx(nnn+1)=inew(idw)
	            nnn=nnn+1
57	            continue
	         enddo

	         kpoin=kpoin+1
	         mpoin(kpoin)=nnn
	         rhs(kpoin)=0.
c==============================================
	      elseif((nsouth.le.i2).and.(i2.le.n))then
c--if i2 is in the s hemisphere (exclude the blocks bounded to the n
c--by the equator because the variation across the equator has already
c--been accounted for).
	         call coordsuper(i2,rila,rilo,nsqrs,nlatzones,eq_incr)
	         upla=rila+eq_incr
c	         ilol=(rilo-rloin_2)*100.+1
c	         ilor=(rilo+rloin_2)*100.-1
	         rilol=rilo-rloin_2
	         rilor=rilo+rloin_2
c	         iupla=upla*100.
c	         iupl=superisqre(iupla,ilol,
c     &nsqrs,nsqtot,nlatzones,n,eq_incr)
c	         iupr=superisqre(iupla,ilor,
c     &nsqrs,nsqtot,nlatzones,n,eq_incr)
	         iupl=isqre(upla,rilol,nsqrs,nsqtot,nlatzones,n,eq_incr)
	         iupr=isqre(upla,rilor,nsqrs,nsqtot,nlatzones,n,eq_incr)
c--define the row corresponding to
c--variation wrt latitude (s hemisphere)
	         g(nnn+1)=weight
	         indx(nnn+1)=inew(i2)
	         nnn=nnn+1

	         peso=1./(iupr-iupl+1)
	         do iup=iupl,iupr
c***
	            do k=1,kireso
	               if(iup.eq.iresolu(k))then
	                  call rang(idw,xlamin,xlamax,xlomin,xlomax,
     &	nsqrs,nsqtot,nlatzones,n,eq_incr)
c--block iup is in the high-resolution area; 
c--find the corresponding southernmost finer blocks:
	                  fin_incr=(360./nsqrs(i1-1))/ifa
c	                  ilauphi=rila*100.+((eq_incr+eq_hi)/2.)*100.
	                  rlauphi=rila+((eq_incr+eq_hi)/2.)
c	                  ilostart=(xlomin+(fin_incr/2.))*100.
			  rlostart=(xlomin+(fin_incr/2.))
	                  do kuks=1,ifa
c	                     ilouphi=ilostart+((kuks-1)*(fin_incr))*100.
	                     rlouphi=rlostart+(kuks-1)*(fin_incr)
	                     iuphi=isqre(rlauphi,rlouphi,nsqrsh,
     &                  nsqtoth,nlatzohi,numhi,eq_hi)
	                     indx(nnn+1)=inewh(iuphi)
	                     g(nnn+1)=-weight*peso/ifa
	                     nnn=nnn+1
	                  enddo
	                  goto67
	               endif
	            enddo
c***
	            g(nnn+1)=-weight*peso
	            indx(nnn+1)=inew(iup)
	            nnn=nnn+1
67	            continue
	         enddo

	         kpoin=kpoin+1
	         mpoin(kpoin)=nnn
	         rhs(kpoin)=0.
c==============================================
c==============================================
	      endif
47	      continue
	   enddo
	enddo
	print*,'number of rows in the g.d.matrix:',kpoin-kpoin0
	return
	end
c*************
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

c*************
	subroutine coordsuper(nbloc,blocla,bloclo,
     &           nsqrs,nlatzones,eq_incr)
c--given a cell index on the earth's surface, finds longitude and latitude
c--(not colatitude) of its center.
c	dimension nsqrs(nlatzones)
	parameter(nlatzomax=180,nlatzhmax=720)
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
	dimension nsqrsh(nlatzhmax),nsqtoth(nlatzhmax+1)

	ntot=0
c--loop(s) over all the blocks
	do 500 ila=1,nlatzones
c--increment latitude
	   rlati=90.-(eq_incr*(ila-1))
c--calculate increment in longitude for this band
	   rinlo=(360./nsqrs(ila))
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
c*************

	subroutine rang(nsq,xlamin,xlamax,xlomin,xlomax,
     &	nsqrs,nsqtot,nlatzones,n,eq_incr)
c	dimension nsqrs(nlatzones),nsqtot(nlatzones)
	parameter(nlatzomax=180,nlatzhmax=720)
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
	dimension nsqrsh(nlatzhmax),nsqtoth(nlatzhmax+1)
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
c*************
	subroutine parametra2(eq_incr,nsqrs,nsqtot,nsqrsh,nsqtoth,
     &	nlatzones,eq_hi,nlatzohi,ifa,iswit,refgrid)
	parameter(nlatzomax=180,nlatzhmax=720)
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
	dimension nsqrsh(nlatzhmax),nsqtoth(nlatzhmax+1)
c	dimension nsqrs(nlatzones),nsqtot(nlatzones+1)
c	dimension nsqrsh(nlatzohi),nsqtoth(nlatzohi+1)
	parameter(pi=3.1415926536)
	numto=0
	numhi=0
	colat=-eq_incr/2.
	do k=1,nlatzones
c--increment colatitude (and therefore latitude) of the node
	   colat=colat+eq_incr
	   theta=(colat/180.)*pi
c--for this latitudinal zone, compute number of blocks (nsqrs)
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
c	      nsqrs(k)=nsqrs(k)-1
	      goto100
	     else
	     endif
	    elseif(360./nsqrs(k).lt.refgrid)then
101	     if(mod(refgrid,360./nsqrs(k)).ne.0)then
c	      nsqrs(k)=nsqrs(k)+1
	      nsqrs(k)=nsqrs(k)-1
	      goto101
	     else
	     endif
	    endif
	   endif
c----------------------------------

c--take care of finer grid:
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
	print*,'numto=',numto
	print*,"total number of blocks:",nsqtot(nlatzones+1)
	return
	end
c*************

	subroutine numbers(icoar,ifine,westbo,eastbo,southbo,rthnobo,
     &	nlatzones,nlatzohi,nsqrs,nsqtot,nsqrsh,nsqtoth,eq_incr,
     &	eq_hi,numto,numhi,ifa,iresolu,kireso,ire2,kire2)

c	dimension nsqrs(nlatzones),nsqtot(nlatzones+1)
c	dimension nsqrsh(nlatzohi),nsqtoth(nlatzohi+1)
	parameter(nlatzomax=180,nlatzhmax=720)
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
	dimension nsqrsh(nlatzhmax),nsqtoth(nlatzhmax+1)
	dimension iresolu(10000),ire2(10000)

c	print*,nsqrs
	print*,"total number of coarse blocks:",nsqtot(nlatzones+1)
c	print*,nsqrsh
	print*,"total number of fine blocks:",nsqtoth(nlatzohi+1)

c--determine indexes of coarse and fine blocks 
c--within the high resolution area:
	kireso=0
	kire2=0
	do parall=rthnobo,southbo,-eq_incr
	   do rmerid=westbo,eastbo,eq_incr
	      kireso=kireso+1
c	      ilat=parall*100.+0.5
c	      ilon=rmerid*100.+0.5
	      iresolu(kireso)=
     &isqre(parall,rmerid,nsqrs,nsqtot,nlatzones,numto,eq_incr)
c     &superisqre(ilat,ilon,nsqrs,nsqtot,nlatzones,numto,eq_incr)
	      icoarse=iresolu(kireso)
c--finer grid
	      call rang(icoarse,xlamin,xlamax,
     &	xlomin,xlomax,nsqrs,nsqtot,nlatzones,n,eq_incr)
	      do ifila=1,ifa
	      do ifilo=1,ifa
	         kire2=kire2+1
	         xlafi=xlamin+((xlamax-xlamin)/ifa)*(ifila-0.5)
	         xlofi=xlomin+((xlomax-xlomin)/ifa)*(ifilo-0.5)
c	         ilat=xlafi*100.+0.5
c	         ilon=xlofi*100.+0.5
	         ire2(kire2)=
c     &	superisqre(ilat,ilon,nsqrsh,nsqtoth,nlatzohi,numhi,eq_hi)
     &	isqre(xlafi,xlofi,nsqrsh,nsqtoth,nlatzohi,numhi,eq_hi)
	      enddo
	      enddo
	   enddo
	enddo
	print*,kireso,' nonzero elements in array iresolu'
	print*,kire2,' nonzero elements in array ire2'
c--count
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
	icoar=numto-icoar
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

	print*,icoar,ifine
	return
	end
c************************************************************
	subroutine correspc(iresolu,kireso,n,ioffset,inew,iold,ico)
	dimension iresolu(10000),inew(50000),iold(50000)
	ico=0
	do i=1,n
	   do k=1,kireso
	   if(i.eq.iresolu(k))then
	      ico=ico+1
	      inew(i)=-1
c	      print*,i,inew(i)
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
	dimension iresolu(10000),inew(50000),iold(50000)
	ico=0
c	write(44,*)'offset=',ioffset
	do i=1,n
	   do k=1,kireso
	   if(i.eq.iresolu(k))then
	      inew(i)=(i+ioffset)-ico
c	      write(44,*)i,inew(i)
	      goto42
	   endif
	   enddo
	   ico=ico+1
	   inew(i)=-1
c	   write(44,*)i,inew(i)
42	   iold(inew(i))=i
	enddo
	return
	end
c************************************************************
	subroutine rough_d_inside(g,indx,mpoin,nblo,weight,
     &  kpoin,nnn,nonz,
     &	m,rhs,nsqrs,nsqtot,nlatzones,eq_incr,n,
     &	nsqrsco,nsqtotco,nlatzoco,inew,inewh,ifa,co_incr)

	dimension indx(nonz),g(nonz),mpoin(0:m),rhs(m)
	parameter(nlatzomax=180,nlatzhmax=720)
	dimension nsqrsco(nlatzomax),nsqtotco(nlatzomax+1)
	dimension nsqrs(nlatzhmax),nsqtot(nlatzhmax+1)
c	dimension nsqrs(nlatzones),nsqtot(nlatzones+1)
c	dimension nsqrsco(nlatzoco),nsqtotco(nlatzoco+1)
	dimension inewh(50000),inew(50000)

	kpoin0=kpoin
	print*,'gradamp=',weight

	nnorth=n/2
c	print*,'nnorth=',nnorth
	indexch=(nlatzones/2)+1
	nsouth=nnorth+nsqrs(indexch)+1
c	print*,'nsouth=',nsouth

c--loop over all the blocks (consists of 2 loops: latitude and longitude)
	do i1=1,nlatzones
c--rloin_2=half the longitudinal size of the block
	   rloin_2=180./nsqrs(i1)
	   ifirst=nsqtot(i1)+1
	   ilast=nsqtot(i1+1)
	   do i2=ifirst,ilast
c***
	    if(inewh(i2).eq.-1)goto64
c***
c--find iright
	      iright=i2+1
	      if(i2.eq.ilast)iright=ifirst
c--define the row of the damping matrix corresponding to 
c--longitude variation.
c***
	      if(inewh(iright).eq.-1)then
	       call coordsuper(iright,rila,rilo,nsqrs,nlatzones,eq_incr)
c	       ila=rila*100.
c	       ilo=rilo*100.
c	       iright=superisqre(ila,ilo,nsqrsco,nsqtotco,nlatzoco,n,co_incr)
	       iright=isqre(rila,rilo,nsqrsco,nsqtotco,
     &              nlatzoco,n,co_incr)
	       indx(nnn+1)=inewh(i2)
	       indx(nnn+2)=inew(iright)
	       g(nnn+1)=weight/ifa
	       g(nnn+2)=-weight/ifa
	       nnn=nnn+2
	      else
c***

	       indx(nnn+1)=inewh(i2)
	       indx(nnn+2)=inewh(iright)
	       g(nnn+1)=weight
	       g(nnn+2)=-weight
	       nnn=nnn+2
c***
	      endif
c***
	      kpoin=kpoin+1
	      mpoin(kpoin)=nnn
	      rhs(kpoin)=0.


c==============================================
c==============================================
c--find iup and idw
	      if((1.le.i2).and.(i2.le.nnorth))then
c--if i2 is in the n hemisphere
	         call coordsuper(i2,rila,rilo,nsqrs,nlatzones,eq_incr)
	         dwla=rila-eq_incr

cTEST
c		 print*,rila,dwla

	         ilol=(rilo-rloin_2)*100.+1 ! will be useless
	         ilor=(rilo+rloin_2)*100.-1 ! will be useless

		 rilol=rilo-rloin_2
		 rilor=rilo+rloin_2

cTEST
c		 print*,ilol,ilor,rilol,rilor

	         idwla=dwla*100.
c	         idwl=superisqre(idwla,ilol,
c     &nsqrs,nsqtot,nlatzones,n,eq_incr)
c	         idwr=superisqre(idwla,ilor,
c     &nsqrs,nsqtot,nlatzones,n,eq_incr)

		 idwl=isqre(dwla,rilol,nsqrs,nsqtot,nlatzones,n,eq_incr)
		 idwr=isqre(dwla,rilor,nsqrs,nsqtot,nlatzones,n,eq_incr)

cTEST
c		 print*,i2,idwla,idwl,idwr
c		 pause

c--define the rows corresponding to
c--variation wrt latitude (n hemisphere)

	         peso=1./(idwr-idwl+1)

c***
	         iflig=0
	         do idw=idwl,idwr
	          if(inewh(idw).eq.-1)then
	           call coordsuper(idw,rila,rilo,nsqrs,nlatzones,eq_incr)
	           ila=rila*100.
	           ilo=rilo*100.
c	           idw1=superisqre(ila,ilo,nsqrsco,
c     &	nsqtotco,nlatzoco,n,co_incr)
	           idw1=isqre(rila,rilo,nsqrsco,
     &	nsqtotco,nlatzoco,n,co_incr)
	           if(iflig.ne.idw1)then

c	           g(nnn+1)=weight/ifa
	           g(nnn+1)=weight
	           indx(nnn+1)=inewh(i2)
	           nnn=nnn+1
	           indx(nnn+1)=inew(idw1)
c	           g(nnn+1)=-weight/ifa
	           g(nnn+1)=-weight
	           nnn=nnn+1
	           kpoin=kpoin+1
	           mpoin(kpoin)=nnn
	           rhs(kpoin)=0.

	           iflig=idw1
	           endif
	          else
	           g(nnn+1)=weight*peso
	           indx(nnn+1)=inewh(i2)
	           nnn=nnn+1
	           g(nnn+1)=-weight*peso
	           indx(nnn+1)=inewh(idw)
	           nnn=nnn+1
	           kpoin=kpoin+1
	           mpoin(kpoin)=nnn
	           rhs(kpoin)=0.
	          endif
c***
	         enddo

c==============================================
	      elseif((nsouth.le.i2).and.(i2.le.n))then
c--if i2 is in the s hemisphere (exclude the blocks bounded to the n
c--by the equator because the variation across the equator has already
c--been accounted for).
	         call coordsuper(i2,rila,rilo,nsqrs,nlatzones,eq_incr)
	         upla=rila+eq_incr
	         ilol=(rilo-rloin_2)*100.+1
	         ilor=(rilo+rloin_2)*100.-1
	         rilol=(rilo-rloin_2)
	         rilor=(rilo+rloin_2)
	         iupla=upla*100.
c	         iupl=superisqre(iupla,ilol,
c     &nsqrs,nsqtot,nlatzones,n,eq_incr)
c	         iupr=superisqre(iupla,ilor,
c     &nsqrs,nsqtot,nlatzones,n,eq_incr)
	         iupl=isqre(upla,rilol,nsqrs,nsqtot,nlatzones,n,eq_incr)
	         iupr=isqre(upla,rilor,nsqrs,nsqtot,nlatzones,n,eq_incr)
c--define the row corresponding to
c--variation wrt latitude (s hemisphere)

	         peso=1./(iupr-iupl+1)
c***
	         iflig=0
	         do iup=iupl,iupr
	          if(inewh(iup).eq.-1)then
	           call coordsuper(iup,rila,rilo,nsqrs,
     &                   nlatzones,eq_incr)
	           ila=rila*100.
	           ilo=rilo*100.
c	           iup1=superisqre(ila,ilo,nsqrsco,
c     &	nsqtotco,nlatzoco,n,co_incr)
	           iup1=isqre(rila,rilo,nsqrsco,
     &                   nsqtotco,nlatzoco,n,co_incr)
	           if(iflig.ne.iup1)then

c	           g(nnn+1)=weight/ifa
	           g(nnn+1)=weight
	           indx(nnn+1)=inew(iup1)
	           nnn=nnn+1
	           indx(nnn+1)=inew(iup1)
c	           g(nnn+1)=-weight/ifa
	           g(nnn+1)=-weight
	           nnn=nnn+1
	           kpoin=kpoin+1
	           mpoin(kpoin)=nnn
	           rhs(kpoin)=0.

	           iflig=iup1
	           endif
	          else
	           g(nnn+1)=weight*peso
	           indx(nnn+1)=inewh(i2)
	           nnn=nnn+1
	           g(nnn+1)=-weight*peso
	           indx(nnn+1)=inewh(iup)
	           nnn=nnn+1
	           kpoin=kpoin+1
	           mpoin(kpoin)=nnn
	           rhs(kpoin)=0.
	          endif
c***

	         enddo
	         kpoin=kpoin+1
	         mpoin(kpoin)=nnn
	         rhs(kpoin)=0.
c==============================================
c==============================================
	      endif
64	      continue
c--end of latitudinal and longitudinal loops:
	   enddo
	enddo

	print*,'number of rows in the g.d.matrix:',kpoin-kpoin0
	return
	end
