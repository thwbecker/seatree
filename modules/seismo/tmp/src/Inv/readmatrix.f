	subroutine readmatrix(ndata,namexxx,nameind,namepoi,namerhs,
     &mpoin,t,indx,values,icol,jj,nonz,n,m,relwei,cutoff,cutoffn,itcount,nmax,ritcount)

c----this version allows for downweighting of outliers

	dimension t(m),indx(nonz),values(nonz),mpoin(0:m)
        dimension itcount(nmax),ritcount(nmax)
	character*80 namexxx,nameind,namepoi,namerhs
cTEST
	print*,'opening files ',namexxx,nameind,namepoi,namerhs
c        pause
c--read in the matrix
	open(1,file=namexxx,status='old',access='direct',recl=4,form='unformatted')
	open(4,file=nameind,status='old',access='direct',recl=4,form='unformatted')
	open(3,file=namepoi,status='old')
	open(77,file=namerhs,status='old')
	print*,'start from',icol,jj
        print*,ndata,"observations to read"
	icol0=icol
	nrec0=1
	do icol=icol0,icol0+ndata-1
	   read(3,*)mptemp
	   mpoin(icol)=mptemp+mpoin(icol0-1)
	   read(77,*)t(icol)
cTEST
c        print*,mptemp,t(icol)
c--assign additional weight to penalize large anomalies (outliers)
           adtime=abs(t(icol))
           wei2=1.
           if(adtime.gt.cutoff)then
              wei2=exp(cutoff-adtime)
           endif
           if(adtime.lt.cutoffn)then
              wei2=exp(-cutoffn+adtime)
           endif
	   t(icol)=t(icol)*relwei*wei2
c	   t(icol)=t(icol)*relwei
c        print*,mpoin(icol0-1),mpoin(icol),t(icol)!test
	   do nrec=nrec0,mptemp
	      read(4,rec=nrec)indx(jj)
              itcount(indx(jj))=itcount(indx(jj))+1
	      read(1,rec=nrec)values(jj)
              values(jj)=values(jj)*relwei
              ritcount(indx(jj))=ritcount(indx(jj))+values(jj)
cTEST
c        print*,nrec,indx(jj),values(jj)
c        pause
	      jj=jj+1
	   enddo
cTEST
c        pause
	   nrec0=mptemp+1
c	   nrec0=mpoin(icol)+1
	   if(mod(icol,5000).eq.0)print*,icol," rows read"!test
	enddo
c	icol=icol-1
c	print*,"out of loop:",icol!test
	print*,'total number of data so far=',icol,"out of maximum",m
	print*,"total number of nonzero entries=",jj-1,"out ot maximum",nonz
	print*,"average number of matrix entries per datum=",float(jj-1)/float(icol)
cTEST
c        print*,"test:"
c        print*,jj,indx(jj-1),indx(jj)
c        print*,jj,values(jj-1),values(jj)
	close(1)
	close(4)
	close(3)
	close(77)
	return
	end
