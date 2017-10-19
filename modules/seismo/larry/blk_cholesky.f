	parameter(ncoef=1988)
c        parameter(ncoef=12406)
        parameter(ndimata=(ncoef*(ncoef+1))/2)
        character fileout*200,answ*1,txtfile*200,wave_type*1
	character filenm*200,filenv*200
        dimension ata(ndimata),atd(ncoef),x(ncoef),ytemp(ncoef)
        dimension gstore(ndimata)
	parameter(pi=3.141592653589)
	parameter(radian=pi/180.0)
        character chn*6
cTEST
	print*,"size of problem=",ndimata,ncoef

c	print*,"root of name of matrix files?"
c	read*,filen
        print*,"data file used to build the matrix?"
	read*,txtfile
	print*,"file with matrix ATA"
	read*,filenm
	print*,"file with vector ATd"
	read*,filenv
	print*,"name of solution model file?"
	read*,fileout
	print*,"norm damp?"
	read*,dampn
	print*,"roughness damp?"
	read*,dampr
	print*,'remove the average from the solution?'
	read*,answ
c-----------------find reference velocity
	open(2,file=txtfile,status='old')
	read(2,*)wave_type,period,velo0
c--qui sopra usa format 1003 se cosi non funziona!
	czero=velo0/(6371.*radian)
	print*,wave_type,period,velo0,czero
	close(2)
1003	format(8x,a1,8x,f7.4,10x,f8.5)
c-----------------open and read ata and atd files
        open(61,file=filenm,access='direct',
     &  recl=4,form='unformatted',status="old")
        open(62,file=filenv,access='direct',
     &  recl=4,form='unformatted',status="old")
        do kata=1,ndimata
           read(61,rec=kata)ata(kata)
        enddo
        do ib=1,ncoef
	   read(62,rec=ib)atd(ib)
	enddo
	close(61)
	close(62)
c-----------------norm damping
        if(dampn.ne.0.)then
	kata=0
	do i=1,ncoef
	do j=1,i
	kata=kata+1
	if(i.eq.j)ata(kata)=ata(kata)+dampn
	enddo
	enddo
        endif
c-----------------roughness damping
        if(dampr.ne.0.)then
           write(chn,'(i6.6)')ncoef
           open(71,file='gtg.'//chn,access='direct',
     &          recl=4,form='unformatted',status="old")
           open(72,file='gtd.'//chn,access='direct',
     &          recl=4,form='unformatted',status="old")
           do kata=1,ndimata
              read(71,rec=kata)gtg
              ata(kata)=ata(kata)+dampr*gtg
           enddo
           do ib=1,ncoef
              read(72,rec=ib)gtd
              atd(ib)=atd(ib)+dampr*gtd
           enddo
           close(71)
           close(72)
        endif
c-----------------cholesky factorization
	print*,"cholesky factorization on one processor..."
        call choles(ata,gstore,atd,ytemp,x,ncoef,ndimata,nono)
	if(nono.ne.0)then
	print*,"factorization impossible"
	print*,nono
	stop
	endif
	print*,"...OK"

c--remove the average
	ave=0.
	if((answ.eq.'Y').or.(answ.eq.'y'))THEN
	   print*,'remove the average'
	   do i=1,ncoef
	      ave=ave+x(i)
	   enddo
	   ave=ave/ncoef
	   print*,'average=',(-czero*ave*100.),czero,ave
	   open(72,file="averages.txt",access="append")
	   write(72,*)period,(-czero*ave*100.),
     &          filenm,dampn,dampr," CHY"
	   close(72)
	else
	   print*,'keep the average'
	endif

        open(2,file=fileout)
        do i=1,ncoef
           result=1./(1./czero+(x(i)-ave))-czero
           result=result/czero
           result=result*100.
           write(2,*)i,result
c           write(2,*)i,-czero*(x(i)-ave)*100.
        enddo
        close(2)

        end

c------------------------single-processor Cholesky factorization
      subroutine choles(a,g,b,y,x,n,nata,nono)
      implicit real*4 (a-h, o-z)
      real*4 a(nata),g(nata)
      real*4 b(n),y(n),x(n)
c        a= row-wise p.d. symm. system  n*(n+1)/2
c        g= cholesky storage
c        b= r.h.s. vector               n
c        y= temp. vector
c        x= answer vector
c        n= system dimension
c        nono .gt. 0 is the level at which p.d. failed
c        (a,g) and (b,y,x) may be equivalenced.
c----------------------------------------------------------
cTEST
	print*,"size of problem=",n,nata
c-----first compute cholesky decomposition
      nono=0      
      if(a(1).le.0.) then
      nono=1
      return
      endif
      
      g(1)=sqrt(a(1))
      y(1)=b(1)/g(1)

      do 400 i=2,n
      
      kz=(i*(i-1))/2
      g(kz+1)=a(kz+1)/g(1)
      sg=g(kz+1)**2
      y(i)=b(i)-g(kz+1)*y(1)
      
      if(i.gt.2) then
      
      jmax=i-1
      
      do 200 j=2,jmax
      
      gkz=a(kz+j)
      kj=(j*(j-1))/2
      kmax=j-1
      
      do 100 k=1,kmax
      gkz=gkz-g(kz+k)*g(kj+k)
100   continue

      g(kz+j)=gkz/g(kj+j)
      y(i)=y(i)-g(kz+j)*y(j)
      sg=sg+g(kz+j)**2
      
200   continue

      endif
      
      gkz=a(kz+i)-sg
      
      if(gkz.le.0.) then
      nono=i
      return
      endif
      
      g(kz+i)=sqrt(gkz)
      y(i)=y(i)/g(kz+i)
      
400   continue

      kz=(n*(n-1))/2
      x(n)=y(n)/g(kz+n)
      if(n.le.1) return

c-----
c     compute solution for particular rhs
      
      do 600 k=2,n
      
      i=n+1-k
      x(i)=y(i)
      jmin=i+1
      
      do 500 j=jmin,n
      kj=(j*(j-1))/2
      x(i)=x(i)-g(kj+i)*x(j)
500   continue

      kz=(i*(i+1))/2
      x(i)=x(i)/g(kz)
      
600   continue

      return
      end
