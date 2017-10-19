	subroutine parametriz(nlay,rboundary,rescale,iswiresc)
c--defines the re-scaling vector and writes it to a file (rescale.d)
	implicit double precision(a-h,o-z)
	real*8 rboundary(0:20)
	real*4 rescale(20)
	common/cmbradius/rcmb
	parameter(pi=3.141592653589)

c--define the model
	yncr=(6371.-dble(rcmb))/nlay
	do iyn=1,nlay
	   rboundary(iyn)=6371.-(yncr*dble(float(iyn)))
	enddo
c--this way the 1st element of rboundary is the uppermost boundary between
c--layers(the first below the surface), while the last element is the core- 
c--mantle boundary.
	rboundary(0)=6371.
c--0th element is Earth's radius.

c--define the vector used to re-scale:
c--(remember to change its length if we change the n. of layers in the
c--model)
	vrnormal=rboundary(nlay)
     &	          +((rboundary(nlay-1)-rboundary(nlay))/2.)
	vrnormal=(vrnormal*pi*5./180.)**2
	vrnormal=vrnormal*(rboundary(nlay-1)-rboundary(nlay))
	do iresc=1,nlay
	   tempor=rboundary(iresc)
     &	          +((rboundary(iresc-1)-rboundary(iresc))/2.)
	   tempor=(tempor*pi*5./180.)**2
	   tempor=tempor*(rboundary(iresc-1)-rboundary(iresc))
	   rescale(nlay-iresc+1)=sqrt(tempor/vrnormal)
c	   rescale(nlay-iresc+1)=tempor/vrnormal
c	   print*,nlay-iresc+1,rescale(nlay-iresc+1)
	enddo
c=========================================================
c--write rescale on file.
	open(54,file='rescale.d')
	do iresc=1,nlay

C--------------turn off rescaling?-----------
	   IF(ISWIRESC.EQ.0)RESCALE(IRESC)=1.
c--------------------------------------------

	   write(54,*)iresc,rescale(iresc)
	enddo
	close(54)

	RETURN
	END
