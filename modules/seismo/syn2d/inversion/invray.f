      character*300 filen,fileout,filesyn
      parameter(nmax=10000,m=20000,nonz=2000000)
      dimension v(nmax),w(nmax),x(nmax),se(nmax),aty(nmax)
      dimension indx(nonz),values(nonz)
      dimension t(m),ttt(m),mpoin(0:m),hitcount(nmax)
      print*,"name of matrix files (without suffix)?"
      read*,filen
      print*,"name of synthetic data file?"
      read*,filesyn
      print*,"number of observations?"
      read*,nobs
      print*,"horizontal increment of grid?"
      read*,x_incr
      print*,"vertical increment of grid?"
      read*,y_incr
      print*,"horizontal length of gridded region?"
      read*,xtot
      print*,"vertical length of gridded region?"
      read*,ytot
      nxtot=xtot/x_incr
      nytot=ytot/y_incr
      print*,"then ",nxtot," horizontal intervals, and "
     &,nytot," vertical intervals"
      n=nxtot*nytot
      print*,n,"pixels"
      print*,"name of output file?"
      read*,fileout
      print*,"damping parameter?"
      read*,damp

      do k=1,300
         if(filen(k:k).eq." ")goto 74
      enddo
 74   kend=k-1

      
      open(1,file=filen(1:kend)//'.xxx',status='old',access='direct',
     &     recl=4,form='unformatted')
      open(4,file=filen(1:kend)//'.ind',status='old',access='direct',
     &     recl=4,form='unformatted')

      open(3,file=filen(1:kend)//'.pnt',status='old')

      open(7,file=filesyn,status="old")
      do k=1,n
	   hitcount(k)=0
	enddo
c--read matrix
      mpoin(0)=0
      do icol=1,nobs
         read(3,*,end=9)mpoin(icol)
         read(7,*,end=7)t(icol)
         if(mod(icol,1000).eq.0)print*,icol,mpoin(icol),t(icol)
         do jj=mpoin(icol-1)+1,mpoin(icol)
            read(4,rec=jj)indx(jj)
            read(1,rec=jj)values(jj)
c            print*,jj,indx(jj),values(jj)!test
            hitcount(indx(jj))=hitcount(indx(jj))+1
            if(indx(jj).gt.n)then
               print*,jj,indx(jj),values(jj)
               stop
            endif
         enddo
c	pause!test

      enddo
 7    print*,"rhs was shorter or m too small"
 9    print*,'number of rows=',icol-1
      kp=icol-1
      kp0=kp
      nz=jj-1
      close(1)
      close(4)
      close(3)
      close(7)
	open(99,file="hitcount.dat")
	do k=1,n
	   write(99,*)k,hitcount(k)
	enddo
	close(99)
c--save r.h.s.; we are going to need it to compute varce reduction.
      print*,'save r.h.s.',kp0
      do it=1,kp0
         ttt(it)=t(it)
      enddo
	print*,"OK"
c---------------define inversion parameters
	atol=0.
        btol=0.
        conlim=0.
        itnlim=nmax
c        itnlim=n
        nout=1

c----------------------------call lsqr routine
      call lsqr( kp0, n, damp,
     &     1, n, iw, aty,
     &     t, v, w, x, se,
     &     ATOL, BTOL, CONLIM, ITNLIM, NOUT,
     &     ISTOP, ITN, ANORM, ACOND, RNORM, ARNORM, XNORM,
     &     indx,values,mpoin,nonz)

c--compute variance reduction
      print*,'compute variance reduction'
      rnumer=0.
      denom=0.
      do l=1,kp0
         denom=denom+(ttt(l)*ttt(l))
         tot=0.
         do ll=mpoin(l-1)+1,mpoin(l)
            tot=tot+(values(ll)*x(indx(ll)))
c            print*,ll,values(ll),indx(ll),x(indx(ll))!test
         enddo
c         write(*,*)l,ttt(l),tot!test
c         pause!test
         rnumer=rnumer+(ttt(l)-tot)**2
      enddo
      varred=1.-(rnumer/denom)
      write(*,*)'variance reduction=',varred
      rms=0.
      open(1,file=fileout)
      do jj=1,n
c         write(1,*)jj,x(jj)-ave
         write(1,*)jj,x(jj)
         rms=rms+x(jj)*x(jj)
      enddo
      rms=rms/n
      rms=sqrt(rms)
c      write(1,*)'iter.s:',itn,' v.red.:',varred,' damp:',
c     &     damp0,gradamp,gradamp_fin
      close(1)
c
c
c     
      do k=1,80
         if(fileout(k:k).eq." ")goto77
      enddo
 77   open(1,file=fileout(1:k-1)//".xyz")
c log gile
      open(84,file=fileout(1:k-1)//".log")

      write(84,*)'iter.s:',itn,' v.red.:',varred,' damp:',
     &     damp0,gradamp,gradamp_fin
      close(84)
      open(66,file="solstat.log")
      write(66,*)rms,varred
      close(66)
      open(2,file="hitcount.xyz")
      jj=0
      do jy=1,nytot
         do jx=1,nxtot
            jj=jj+1
            write(1,*)
     &(jx-1)*x_incr+x_incr/2.,(jy-1)*y_incr+y_incr/2.,x(jj)
            write(2,*)
     &(jx-1)*x_incr+x_incr/2.,(jy-1)*y_incr+y_incr/2.,hitcount(jj)
	    if(inpx(jx,jy,nxtot,nytot).ne.jj)then
	       print*,inpx(jx,jy,nxtot,nytot),jj," there is a problem"
	    endif
         enddo
      enddo
      close(1)
      close(2)
c---------------rms and misfit
      open(1,file="lcurve.txt",access="append")
      write(1,*)rms,1.-varred
      end
c==================================================================
	subroutine aprod(mode,m,n,x,y,lenin,lenva,iw,aty,
     &       indx,values,mpoin,nonz)
c	implicit real*8(a-h,o-z)
	integer iw(lenin)
	real*4 x(n),y(m),aty(n)
	dimension indx(nonz),values(nonz),mpoin(0:m)

	if(mode.eq.1)then
cTEST
c	print*,"in aprod 1",m
c--compute y=y+A*x
	do k=1,m
	   pro=0.
	   do j=mpoin(k-1)+1,mpoin(k)
cTEST
c	print*,k,mpoin(k-1),mpoin(k),j,values(j),indx(j)
	      pro=pro+values(j)*x(indx(j))
	   enddo
	   y(k)=y(k)+pro
	enddo
cTEST
c	print*,"leave aprod 1"
	return

	elseif(mode.eq.2)then
cTEST
c	print*,"in aprod 2"
	do i=1,n
	   aty(i)=0.
	enddo
c--compute x=x+(A^t)*y
	do k=1,m
	   do j=mpoin(k-1)+1,mpoin(k)
	      aty(indx(j))=aty(indx(j))+values(j)*y(k)
	   enddo
	enddo
	do i=1,n
	   x(i)=x(i)+aty(i)
	enddo
cTEST
c	print*,"leave aprod 2"
	return

	else
	print*,'error: mode=',mode
	stop
	endif
	end
C================================================================
* From rfischer@seismology.harvard.edu Wed Jul 23 15:30:22 1997
* Retrieved by Bob Fischer from: http://www.netlib.org/linalg/lsqr
* From arpa!sol-michael.stanford.edu!mike 5 May 89 23:53:00 PDT
      SUBROUTINE LSQR  ( M, N, DAMP,
     $                   LENIW, LENRW, IW, rw,
     $                   U, V, W, X, SE,
     $                   ATOL, BTOL, CONLIM, ITNLIM, NOUT,
     $   ISTOP, ITN, ANORM, ACOND, RNORM, ARNORM, XNORM,
     &       indx,values,mpoin,nonz)

c--I had to add next line to pass the matrix to aprod without using a common
	dimension indx(nonz),values(nonz),mpoin(0:m)

c      EXTERNAL           APROD
      INTEGER            M, N, LENIW, LENRW, ITNLIM, NOUT, ISTOP, ITN
      INTEGER            IW(LENIW)
c      DOUBLE PRECISION   RW(LENRW), U(M), V(N), W(N), X(N), SE(N),
      real*4  RW(LENRW), U(M), V(N), W(N), X(N), SE(N),
     $                   ATOL, BTOL, CONLIM, DAMP,
     $                   ANORM, ACOND, RNORM, ARNORM, XNORM

*-----------------------------------------------------------------------
*     Intrinsics and local variables

      INTRINSIC          ABS, MOD, SQRT
      INTEGER            I, NCONV, NSTOP
c      DOUBLE PRECISION   DNRM2
	real*4 dnrm2
c      DOUBLE PRECISION   ALFA, BBNORM, BETA, BNORM,
      real*4   ALFA, BBNORM, BETA, BNORM,
     $                   CS, CS1, CS2, CTOL, DAMPSQ, DDNORM, DELTA,
     $                   GAMMA, GAMBAR, PHI, PHIBAR, PSI,
     $                   RES1, RES2, RHO, RHOBAR, RHBAR1, RHBAR2,
     $                   RHS, RTOL, SN, SN1, SN2,
     $                   T, TAU, TEST1, TEST2, TEST3,
     $                   THETA, T1, T2, T3, XXNORM, Z, ZBAR

c      DOUBLE PRECISION   ZERO,           ONE
      PARAMETER        ( ZERO = 0.,  ONE = 1. )

      CHARACTER*16       ENTER, EXIT
      CHARACTER*60       MSG(0:7)

      DATA               ENTER /' Enter LSQR.    '/,
     $                   EXIT  /' Exit  LSQR.    '/

      DATA               MSG
     $ / 'The exact solution is  X = 0',
     $   'Ax - b is small enough, given ATOL, BTOL',
     $   'The least-squares solution is good enough, given ATOL',
     $   'The estimate of cond(Abar) has exceeded CONLIM',
     $   'Ax - b is small enough for this machine',
     $   'The least-squares solution is good enough for this machine',
     $   'Cond(Abar) seems to be too large for this machine',
     $   'The iteration limit has been reached' /
*-----------------------------------------------------------------------


*     Initialize.

      IF (NOUT .GT. 0)
     $   WRITE(NOUT, 1000) ENTER, M, N, DAMP, ATOL, CONLIM, BTOL, ITNLIM
      ITN    =   0
      ISTOP  =   0
      NSTOP  =   0
      CTOL   =   ZERO
      IF (CONLIM .GT. ZERO) CTOL = ONE / CONLIM
      ANORM  =   ZERO
      ACOND  =   ZERO
      BBNORM =   ZERO
      DAMPSQ =   DAMP**2
      DDNORM =   ZERO
      RES2   =   ZERO
      XNORM  =   ZERO
      XXNORM =   ZERO
      CS2    = - ONE
      SN2    =   ZERO
      Z      =   ZERO

      DO 10  I = 1, N
         V(I)  =  ZERO
         X(I)  =  ZERO
        SE(I)  =  ZERO
   10 CONTINUE

*     Set up the first vectors U and V for the bidiagonalization.
*     These satisfy  BETA*U = b,  ALFA*V = A(transpose)*U.

      ALFA   =   ZERO
      BETA   =   DNRM2 ( M, U, 1 )

      IF (BETA .GT. ZERO) THEN
         CALL DSCAL ( M, (ONE / BETA), U, 1 )
         CALL APROD ( 2, M, N, V, U, LENIW, LENRW, IW, RW ,
     &       indx,values,mpoin,nonz)
         ALFA   =   DNRM2 ( N, V, 1 )
      END IF

      IF (ALFA .GT. ZERO) THEN
         CALL DSCAL ( N, (ONE / ALFA), V, 1 )
         CALL DCOPY ( N, V, 1, W, 1 )
      END IF

      ARNORM =   ALFA * BETA
      IF (ARNORM .EQ. ZERO) GO TO 800

      RHOBAR =   ALFA
      PHIBAR =   BETA
      BNORM  =   BETA
      RNORM  =   BETA

      IF (NOUT   .GT.  0  ) THEN
         IF (DAMPSQ .EQ. ZERO) THEN
             WRITE(NOUT, 1200)
         ELSE
             WRITE(NOUT, 1300)
         END IF
         TEST1  = ONE
         TEST2  = ALFA / BETA
         WRITE(NOUT, 1500) ITN, X(1), RNORM, TEST1, TEST2
         WRITE(NOUT, 1600)
      END IF

*     ------------------------------------------------------------------
*     Main iteration loop.
*     ------------------------------------------------------------------
  100 ITN    = ITN + 1
	print*,'iteration:',itn

*     Perform the next step of the bidiagonalization to obtain the
*     next  BETA, U, ALFA, V.  These satisfy the relations
*                BETA*U  =  A*V  -  ALFA*U,
*                ALFA*V  =  A(transpose)*U  -  BETA*V.

      CALL DSCAL ( M, (- ALFA), U, 1 )
      CALL APROD ( 1, M, N, V, U, LENIW, LENRW, IW, RW,
     &       indx,values,mpoin,nonz)

      BETA   =   DNRM2 ( M, U, 1 )
      BBNORM =   BBNORM  +  ALFA**2  +  BETA**2  +  DAMPSQ

      IF (BETA .GT. ZERO) THEN
         CALL DSCAL ( M, (ONE / BETA), U, 1 )
         CALL DSCAL ( N, (- BETA), V, 1 )
         CALL APROD ( 2, M, N, V, U, LENIW, LENRW, IW, RW ,
     &       indx,values,mpoin,nonz)
         ALFA   =   DNRM2 ( N, V, 1 )
         IF (ALFA .GT. ZERO) THEN
            CALL DSCAL ( N, (ONE / ALFA), V, 1 )
         END IF
      END IF

*     Use a plane rotation to eliminate the damping parameter.
*     This alters the diagonal (RHOBAR) of the lower-bidiagonal matrix.

      RHBAR2 = RHOBAR**2  +  DAMPSQ
      RHBAR1 = SQRT( RHBAR2 )
      CS1    = RHOBAR / RHBAR1
      SN1    = DAMP   / RHBAR1
      PSI    = SN1 * PHIBAR
      PHIBAR = CS1 * PHIBAR

*     Use a plane rotation to eliminate the subdiagonal element (BETA)
*     of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.

      RHO    =   SQRT( RHBAR2  +  BETA**2 )
      CS     =   RHBAR1 / RHO
      SN     =   BETA   / RHO
      THETA  =   SN * ALFA
      RHOBAR = - CS * ALFA
      PHI    =   CS * PHIBAR
      PHIBAR =   SN * PHIBAR
      TAU    =   SN * PHI

*     Update  X, W  and the standard error estimates.

      T1     =   PHI   / RHO
      T2     = - THETA / RHO
      T3     =   ONE   / RHO

      DO 200  I =  1, N
         T      =  W(I)
         X(I)   =  T1*T  +  X(I)
         W(I)   =  T2*T  +  V(I)
         T      = (T3*T)**2
         SE(I)  =  T     +  SE(I)
         DDNORM =  T     +  DDNORM
  200 CONTINUE

*     Use a plane rotation on the right to eliminate the
*     super-diagonal element (THETA) of the upper-bidiagonal matrix.
*     Then use the result to estimate  norm(X).

      DELTA  =   SN2 * RHO
      GAMBAR = - CS2 * RHO
      RHS    =   PHI    - DELTA * Z
      ZBAR   =   RHS    / GAMBAR
      XNORM  =   SQRT( XXNORM    + ZBAR **2 )
      GAMMA  =   SQRT( GAMBAR**2 + THETA**2 )
      CS2    =   GAMBAR / GAMMA
      SN2    =   THETA  / GAMMA
      Z      =   RHS    / GAMMA
      XXNORM =   XXNORM + Z**2

*     Test for convergence.
*     First, estimate the norm and condition of the matrix  Abar,
*     and the norms of  rbar  and  Abar(transpose)*rbar.

      ANORM  =   SQRT( BBNORM )
      ACOND  =   ANORM * SQRT( DDNORM )
      RES1   =   PHIBAR**2
      RES2   =   RES2  +  PSI**2
      RNORM  =   SQRT( RES1 + RES2 )
      ARNORM =   ALFA  * ABS( TAU )

*     Now use these norms to estimate certain other quantities,
*     some of which will be small near a solution.

      TEST1  =   RNORM /  BNORM
      TEST2  =   ZERO
      IF (RNORM .GT. ZERO) TEST2 = ARNORM / (ANORM * RNORM)
      TEST3  =   ONE   /  ACOND
      T1     =   TEST1 / (ONE  +  ANORM * XNORM / BNORM)
      RTOL   =   BTOL  +  ATOL *  ANORM * XNORM / BNORM

*     The following tests guard against extremely small values of
*     ATOL, BTOL  or  CTOL.  (The user may have set any or all of
*     the parameters  ATOL, BTOL, CONLIM  to zero.)
*     The effect is equivalent to the normal tests using
*     ATOL = RELPR,  BTOL = RELPR,  CONLIM = 1/RELPR.

      T3     =   ONE + TEST3
      T2     =   ONE + TEST2
      T1     =   ONE + T1
      IF (ITN .GE. ITNLIM) ISTOP = 7
      IF (T3  .LE. ONE   ) ISTOP = 6
      IF (T2  .LE. ONE   ) ISTOP = 5
      IF (T1  .LE. ONE   ) ISTOP = 4

*     Allow for tolerances set by the user.

      IF (TEST3 .LE. CTOL) ISTOP = 3
      IF (TEST2 .LE. ATOL) ISTOP = 2
      IF (TEST1 .LE. RTOL) ISTOP = 1
*     ==================================================================

*     See if it is time to print something.

      IF (NOUT  .LE.  0       ) GO TO 600
      IF (N     .LE. 40       ) GO TO 400
      IF (ITN   .LE. 10       ) GO TO 400
      IF (ITN   .GE. ITNLIM-10) GO TO 400
      IF (MOD(ITN,10) .EQ. 0  ) GO TO 400
      IF (TEST3 .LE.  2.0*CTOL) GO TO 400
      IF (TEST2 .LE. 10.0*ATOL) GO TO 400
      IF (TEST1 .LE. 10.0*RTOL) GO TO 400
      IF (ISTOP .NE.  0       ) GO TO 400
      GO TO 600

*     Print a line for this iteration.

  400 WRITE(NOUT, 1500) ITN, X(1), RNORM, TEST1, TEST2, ANORM, ACOND
      IF (MOD(ITN,10) .EQ. 0) WRITE(NOUT, 1600)
*     ==================================================================

*     Stop if appropriate.
*     The convergence criteria are required to be met on  NCONV
*     consecutive iterations, where  NCONV  is set below.
*     Suggested value:  NCONV = 1, 2  or  3.

  600 IF (ISTOP .EQ. 0) NSTOP = 0
      IF (ISTOP .EQ. 0) GO TO 100
      NCONV  =   1
      NSTOP  =   NSTOP + 1
      IF (NSTOP .LT. NCONV  .AND.  ITN .LT. ITNLIM) ISTOP = 0
      IF (ISTOP .EQ. 0) GO TO 100
*     ------------------------------------------------------------------
*     End of iteration loop.
*     ------------------------------------------------------------------


*     Finish off the standard error estimates.

      T    =   ONE
      IF (M      .GT.   N )  T = M - N
      IF (DAMPSQ .GT. ZERO)  T = M
      T    =   RNORM / SQRT( T )

      DO 700  I = 1, N
         SE(I)  = T * SQRT( SE(I) )
  700 CONTINUE

*     Print the stopping condition.

  800 IF (NOUT .GT. 0) THEN
         WRITE(NOUT, 2000) EXIT, ISTOP, ITN,
     $                     EXIT, ANORM, ACOND,
     $                     EXIT, RNORM, ARNORM,
     $                     EXIT, BNORM, XNORM
         WRITE(NOUT, 3000) EXIT, MSG(ISTOP)
      END IF

  900 RETURN

*     ------------------------------------------------------------------
 1000 FORMAT(// 1P, A, '  Least-squares solution of  A*x = b'
     $    / ' The matrix  A  has', I7, ' rows   and', I7, ' columns'
     $    / ' The damping parameter is         DAMP   =', E10.2
     $    / ' ATOL   =', E10.2, 15X,        'CONLIM =', E10.2
     $    / ' BTOL   =', E10.2, 15X,        'ITNLIM =', I10)
 1200 FORMAT(// '   Itn       x(1)           Function',
     $   '     Compatible   LS        Norm A    Cond A' /)
 1300 FORMAT(// '   Itn       x(1)           Function',
     $   '     Compatible   LS     Norm Abar Cond Abar' /)
 1500 FORMAT(1P, I6, 2E17.9, 4E10.2)
 1600 FORMAT(1X)
 2000 FORMAT(/ 1P, A, 6X, 'ISTOP =', I3,   16X, 'ITN    =', I9
     $       /     A, 6X, 'ANORM =', E13.5, 6X, 'ACOND  =', E13.5
     $       /     A, 6X, 'RNORM =', E13.5, 6X, 'ARNORM =', E13.5,
     $       /     A, 6X, 'BNORM =', E13.5, 6X, 'XNORM  =', E13.5)
 3000 FORMAT( A, 6X, A )
c     ------------------------------------------------------------------
c     End of LSQR
      end


C======================================================================
C======================================================================

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

C================================================================
C================================================================

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

C===============================================================
C===============================================================

c      DOUBLE PRECISION FUNCTION DNRM2 ( N, X, INCX )
      real FUNCTION DNRM2 ( N, X, INCX )
*     .. Scalar Arguments ..
      INTEGER                           INCX, N
*     .. Array Arguments ..
c      DOUBLE PRECISION                  X( * )
      real*4               X( * )
*     ..
*
*  DNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DNRM2 := sqrt( x'*x )
*
*
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to DLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*
*     .. Parameters ..
c      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1., ZERO = 0. )
*     .. Local Scalars ..
      INTEGER               IX
c      DOUBLE PRECISION      ABSXI, NORM, SCALE, SSQ
      real*4     ABSXI, NORM, SCALE, SSQ
*     .. Intrinsic Functions ..
      INTRINSIC             ABS, SQRT
*     ..
*     .. Executable Statements ..
      IF( N.LT.1 .OR. INCX.LT.1 )THEN
         NORM  = ZERO
      ELSE IF( N.EQ.1 )THEN
         NORM  = ABS( X( 1 ) )
      ELSE
         SCALE = ZERO
         SSQ   = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
*
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( X( IX ).NE.ZERO )THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI )THEN
                  SSQ   = ONE   + SSQ*( SCALE/ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SSQ   = SSQ   +     ( ABSXI/SCALE )**2
               END IF
            END IF
   10    CONTINUE
         NORM  = SCALE * SQRT( SSQ )
      END IF
*
      DNRM2 = NORM
      RETURN
*
*     End of DNRM2.
*
      END
           function inpx(nx,ny,nxtot,nytot)
           inpx=(ny-1)*nxtot+nx
           return
           end

