c
c
c
c
c     shoot rays through a selection of sources.txt and 
c 
c     and record them at all recivers listed in receivers.txt
c
c
c
c this version reads in event-receiver geometry from file
c
      parameter(nmax=500000,isnmax=5000)
      dimension x(nmax),y(nmax),el(nmax),ds(nmax),indnz(nmax)
      dimension source(isnmax,2)
      character filen*80
      
      integer mode,nsource
      
      nsource=0
      nreco=0
      
      print*,"horizontal increment of grid?"
      read*,x_incr
      print*,"vertical increment of grid?"
      read*,y_incr
      print*,"horizontal length of gridded region?"
      read*,xtot
      print*,"vertical length of gridded region?"
      read*,ytot
      print*,"name for output files (suffix will be added)?"
      read*,filen
      do k=1,80
         if(filen(k:k).eq." ")goto 74
      enddo
 74   kend=k-1
 22   print*,"store matrix in binary(1) or ascii(2) format?"
      read*,ifo
      if(ifo.eq.1)then
         open(11,file=filen(1:kend)//'.xxx',access='direct',recl=4,form='unformatted')
         open(22,file=filen(1:kend)//'.ind',access='direct',recl=4,form='unformatted')
      elseif(ifo.eq.2)then
         open(11,file=filen(1:kend)//'.xxx')
         open(22,file=filen(1:kend)//'.ind')
      else
         goto22
      endif

      open(33,file=filen(1:kend)//'.pnt')
      print*,"increment along ray path?"
      read *,d_incr



c--------------------------------------------2d parameterization
      nxtot=xtot/x_incr
      nytot=ytot/y_incr
      print*,"then ",nxtot," horizontal intervals, and "
     &     ,nytot," vertical intervals"
      print*,nxtot*nytot,"pixels"
      if(nxtot*nytot.gt.nmax)stop "too many"

c-------------------------------------main loop starts here
      open(29,file="paths.txt")

      xcen=xtot/2.
      ycen=ytot/2.
      

      
c
c start of main loop until krow == ndata 
c
c

      ysran_old = 1e20
      xsran_old = 1e20
      
      nsource=0
      krow=0

      nrecord=0      
 1    continue
c-------------------------generate a source, a receiver, and the straight path between them


      read(29,*,end=5)xsran,ysran,xrran,yrran ! receivers
      
      if((abs(xsran-xsran_old).gt.1e-6).and.
     &     (abs(ysran-ysran_old).gt.1e-6))then
c     new source
         nsource=nsource+1
c         if(nsource.gt.1)then
c            print *,nrecord,' records for source ',nsource-1
c         endif
         nrecord=1
         xsran_old = xsran
         ysran_old = ysran
      else
         nrecord=nrecord+1
      endif
c     data counter
      krow=krow+1
      
      epidist=sqrt((xsran-xrran)**2+(ysran-yrran)**2)
      
      azi=asin(abs(xsran-xrran)/epidist)
      dx=d_incr*sin(azi)
      dy=d_incr*cos(azi)
      
      if(xrran.lt.xsran)dx=-dx
      if(yrran.lt.ysran)dy=-dy
      
      d=0.
      k=1
      x(1)=xsran
      el(1)=ysran
 3    k=k+1
      x(k)=x(k-1)+dx
      el(k)=el(k-1)+dy
      d=d+d_incr
      if(d.gt.epidist)goto4
      goto3
 4    k=k-1
      x(k)=xrran
      el(k)=yrran
      np=k
      
c------------------------------------------------------project path onto grid:
c------------------------------------------------------identify sampled pixels
c------------------------------------------------------and increment matrix
      x0=x(1)
      y0=el(1)
      xint0=x(1)
      yint0=el(1)
      nx0=(x0/x_incr)+1
      ny0=(y0/y_incr)+1
      ipx0=inpx(nx0,ny0,nxtot,nytot)
      index=0
      
      
      
c------------loop over all points on the ray path
      do k=2,np
         nx=(x(k)/x_incr)+1
         ny=(el(k)/y_incr)+1
         
         if(inpx(nx,ny,nxtot,nytot).ne.ipx0)then ! crossed over to another pixel
            if(ny.ne.ny0)then   !vertical intersection
               yint=min(ny,ny0)*y_incr
               if(ny.lt.ny0)yint=ny*y_incr !DO I NEED THIS LINE?
               a=(x(k)-x0)/(el(k)-y0)
               b=x0-a*y0
               xint=a*yint+b    !find x of intersection by linear interpolation
            elseif(nx.ne.nx0)then !horizontal intersection
               xint=min(nx,nx0)*x_incr
               a=(el(k)-y0)/(x(k)-x0)
               b=y0-a*x0
               yint=a*xint+b    !find y of intersection by linear interpolation
            endif
            index=index+1
            ds(index)=sqrt(((xint-xint0)**2)+((yint-yint0)**2))
c     write(44,*)xint,yint,xint0,yint0,ipx0!test
            indnz(index)=ipx0
            xint0=xint
            yint0=yint
            nx0=nx
            ny0=ny
            ipx0=inpx(nx0,ny0,nxtot,nytot)
         endif
         x0=x(k)
         y0=el(k)
      enddo                     ! ray path is done
c----------------------------------pixel where receiver is
      if(x(np).ne.xint0)then
         index=index+1
         ds(index)=sqrt(((x(np)-xint0)**2)+((el(np)-yint0)**2))
         nx=(x(np)/x_incr)+1
         ny=(el(np)/y_incr)+1           
         indnz(index)=inpx(nx,ny,nxtot,nytot)
      endif
c----------------------------------store sparse matrix row on files
      if(ifo.eq.1)then
         do iw=1,index
            nreco=nreco+1
            write(11,rec=nreco)ds(iw)
            write(22,rec=nreco)indnz(iw)
         enddo
      else
         do iw=1,index
            nreco=nreco+1
            write(11,*)ds(iw)
            write(22,*)indnz(iw)
         enddo
      endif
      write(33,*)nreco
      
      if(mod(krow,500).eq.0)print*,krow
      
      
c-------------------------------------main loop ends here
      goto 1

 5    continue

c      print *,nrecord,' records for source ',nsource-1
      
      print *,'final data number ',krow
      close(11)
      close(22)
      close(29)
      close(33)
      
      end
      
      function inpx(nx,ny,nxtot,nytot)
      inpx=(ny-1)*nxtot+nx
      return
      end
