c--multiply matrix by model vector to obtain synthetic data.
c--tested using a constant model vector =1 so data are lengths of rays.
c--random noise is added to predicted datum
	character*80 matfile,indfile,pntfile,modfile,nparfile,nobsfile
	character*80 vecfile
        parameter(nonz=6000000)
        parameter(m=100000)
	parameter(nmax=30000)
	dimension x(nmax)
        dimension indx(nonz),values(nonz),mpoin(0:m)
        dimension itcount(nmax)
        dimension ritcount(nmax)
	character chnfree*6,choutfile*80
	logical use_random

        print*,"matrix file?"
        read(*,"(a80)")matfile
        print*,"index file?"
        read(*,"(a80)")indfile
        print*,"pointer file?"
        read(*,"(a80)")pntfile
        print*,"model file?"
        read(*,"(a80)")modfile
	print*,"output synthetic data file?"
	read*,choutfile
	print*,"number of free parameters?"
	read*,nfree
        print*,"number of observations?"
        read*,ndata
        print*,"uncertainty for generated data?"
        read*,sigma
	if(abs(sigma) .lt. 1e-7)then
	   use_random=.false.
	else
	   use_random=.true.
	endif
	print*,"seed for random number generation?"
	read*,iseed
c-----------------initialize hitcount arrays (not really needed, but)
        do ihit=1,nmax
           itcount(ihit)=0
           ritcount(ihit)=0
        enddo
	write(chnfree,"(i6.6)")nfree
c-----open and read model file
	open(1,file=modfile,status="old")
	do k=1,nfree
	   read(1,*)kin,x(k)
	   if(kin.ne.k)stop "coefficient index mismatch"
	enddo
c	write(*,*)(x(k),k=1,nfree)
c	pause

	close(1)
c-----open and read A matrix files
	print*,"reading the matrix and the data..."
	relwei=1.
	icol=1
	jj=1
	mpoin(0)=0
        call readmatrix(ndata,matfile,indfile,pntfile,vecfile,
     &mpoin,indx,itcount,values,icol,jj,nonz,nfree,m,relwei,ritcount,
     &nmax)
	open(99,file=choutfile)
        do l=1,icol-1
           tot=0.
           do ll=mpoin(l-1)+1,mpoin(l)
              tot=tot+(values(ll)*x(indx(ll)) )
           enddo
	   if(use_random)then
c---now add random noise to predicted datum. Box-Muller transformation.
	      tot = tot + gasdev(iseed)*sigma
	   endif
	   write(99,*)tot
        enddo
	close(99)
	end

        subroutine readmatrix(ndata,NAMEXXX,NAMEIND,NAMEPOI,NAMERHS,
     &MPOIN,INDX,ITCOUNT,VALUES,ICOL,JJ,NONZ,N,M,relwei,ritcount,nmax)
        dimension itcount(nmax),ritcount(nmax)
        dimension indx(nonz),values(nonz),mpoin(0:m)
        CHARACTER*80 NAMEXXX,NAMEPOI,NAMEIND,NAMERHS
        print*,'opening files ',NAMEXXX,NAMEIND,NAMEPOI,NAMERHS
c--read in the matrix
        open(1,file=namexxx,
     &status='old',access='direct',recl=4,form='unformatted')
        open(4,file=nameind,
     &status='old',access='direct',recl=4,form='unformatted')
        open(3,file=namepoi,status='old')

        print*,'start from',icol,jj
        print*,'weight=',relwei
        icol0=icol-1
        nrec0=1
        do icol=icol0+1,icol0+ndata
        read(3,*,err=153)mptemp
        mpoin(icol)=mptemp+mpoin(icol0)
        do nrec=nrec0,mpoin(icol)
           read(4,rec=nrec)indx(jj)
	   if(indx(jj).gt.n)then
	      print *,"undefined voxel index" !test
	      print *,indx(jj),n
	      stop
	   endif
           read(1,rec=nrec)values(jj)
           itcount(indx(jj))=itcount(indx(jj))+1
           values(jj)=values(jj)*relwei
c	   write(*,*)jj,indx(jj),values(jj),icol,relwei!test
           jj=jj+1
        enddo
c	pause!test

        nrec0=mpoin(icol)+1
        if(mod(icol,1000).eq.0)print*,icol," rows read"
        enddo
        print*,'total number of data so far=',icol-1
        print*,"total number of nonzero entries=",jj-1
        close(1)
        close(4)
        close(3)
        close(77)
        return
153     print*,"error while reading pointer", icol,icol0,ndata,mptemp
        stop
154     print*,"error while reading data vector"
        end



