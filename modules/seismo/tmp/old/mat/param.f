        subroutine param(eq_incr,nsqrs,nsqtot,nlatzones,numto,iswit,refgrid
     &	,nlatzomax)
c---find vectors nsqrs and nsqtot that define a block parameterization
c---eq_incr,iswit,refgrid are input, the rest is output
c        dimension nsqrs(nlatzones),nsqtot(nlatzones+1)
        dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
	parameter(pi=3.1415926536)
	numto=0
	colat=-eq_incr/2.
	do k=1,nlatzones
	   colat=colat+eq_incr
	   theta=(colat/180.)*pi
c--for this latitudinal zone, compute number of blocks (nsqrs)
	   deltalon=eq_incr/(sin(theta))
	   nsqrs(k)=(360./deltalon)+1
           if(mod(nsqrs(K),2).ne.0)nsqrs(k)=nsqrs(k)-1
c--if requested, correct nsqrs(k) so the grid is compatible to reference grid
	   if(iswit.eq.1)then
              if(360./nsqrs(k).ge.refgrid)then
 100             if(mod(360./nsqrs(k),refgrid).ne.0)then
                    nsqrs(k)=nsqrs(k)+1
                    goto100
                 else
                 endif
              elseif(360./nsqrs(k).lt.refgrid)then
 101             if(mod(refgrid,360./nsqrs(k)).ne.0)then
                    nsqrs(k)=nsqrs(k)-1
                    goto101
                 else
                 endif
              endif
	   endif
           if(MOD(NSQRS(K),2).ne.0)stop "nsqrs has to be even"
c---------------------------
	   nsqtot(k)=numto
	   numto=numto+nsqrs(k)
	enddo
	nsqtot(nlatzones+1)=numto
cTEST	
c	do kt=1,nlatzones
c	print*,kt,nsqtot(kt),nsqrs(kt)
c	enddo
	print*,nlatzones+1,nsqtot(nlatzones+1)
	open(29,file="num.px")
	write(29,*)nsqtot(nlatzones+1)
	close(29)
	print*,"total number of blocks:",nsqtot(nlatzones+1)
	return
	end
