        character*50 namexxx,nameind,namepoi,namerhs
        PARAMETER(NONZ=60000000)
        PARAMETER(M=800000)
	parameter(n=4592*4*15)
	dimension indx(nonz),values(nonz),mpoin(0:m)
	integer itcount(n)
	real*4 ritcount(n)
        
	print*,n
	pause

	print*,"a matrix?"
	read(*,*)namexxx
	print*,"index array?"
	read*,nameind
	print*,"pointer array?"
	read*,namepoi
	print*,"data vector?"
	read*,namerhs
	print*,"weight?"
	read*,relwei
	print*,"number of observations in this subset?"
	read*,ndata
	icol=1
	jj=1
	MPOIN(0)=0
	call readmatrix(ndata,NAMEXXX,NAMEIND,NAMEPOI,NAMERHS,
     &MPOIN,T,INDX,ITCOUNT,VALUES,ICOL,JJ,NONZ,N,M,relwei,ritcount)

        end

	subroutine readmatrix(ndata,NAMEXXX,NAMEIND,NAMEPOI,NAMERHS,
     &MPOIN,T,INDX,ITCOUNT,VALUES,ICOL,JJ,NONZ,N,M,relwei,ritcount)

	dimension t(m),itcount(n),ritcount(n)
	dimension indx(nonz),values(nonz),mpoin(0:m)

	CHARACTER*50 NAMEXXX,NAMEPOI,NAMEIND,NAMERHS

	print*,'opening files ',NAMEXXX,NAMEIND,NAMEPOI,NAMERHS

c--read in the matrix
	OPEN(1,FILE=namexxx,
     &status='old',access='direct',recl=4,form='unformatted')
	OPEN(4,FILE=nameind,
     &status='old',access='direct',recl=4,form='unformatted')
	OPEN(3,FILE=namepoi,status='old')
	open(77,file=namerhs,status='old')

	print*,'start from',icol,jj
	print*,'weight=',relwei

	icol0=icol-1
	nrec0=1
	do icol=icol0+1,icol0+ndata
	read(3,*,err=153)mptemp
	mpoin(icol)=mptemp+mpoin(icol0)
	read(77,*,err=154)t(icol)
c--assign weight to rhs:
	t(icol)=t(icol)*relwei
cTEST
	print*,nrec0,mpoin(icol),t(icol),relwei

	do nrec=nrec0,mpoin(icol)
	   read(4,rec=nrec)indx(jj)
	   read(1,rec=nrec)values(jj)
           itcount(indx(jj))=itcount(indx(jj))+1
           values(jj)=values(jj)*relwei
cTEST
	print*,nrec,indx(jj),values(jj),jj
           jj=jj+1
cTEST
	   if(mod(jj,100).eq.0)print*,jj," rows read"
	enddo
	nrec0=mpoin(icol)+1
cTEST
        pause

	enddo
	print*,'total number of data so far=',icol
	print*,"total number of nonzero entries=",jj-1
	close(1)
	close(4)
	close(3)
	close(77)
	return
153	print*,"error while reading pointer", icol,icol0,ndata,mptemp
	stop
154	print*,"error while reading data vector"
	END
