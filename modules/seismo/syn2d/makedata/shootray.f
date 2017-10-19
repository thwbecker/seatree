c
c--generates random distribution of sources and stations. 
c--assuming rays are straight lines, computes corresp. tomographic matrix
c

        parameter(nmax=500000,isnmax=5000)
	dimension x(nmax),y(nmax),el(nmax),ds(nmax),indnz(nmax)
        dimension source(isnmax,2)
        character filen*80
        character chrow*6 !TEST
        real ran2,rtmp
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
        print*,"how many data to build?"
        read*,ndata
	do k=1,80
	   if(filen(k:k).eq." ")goto 74
	enddo
 74	kend=k-1
 22     print*,"store matrix in binary(1) or ascii(2) format?"
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
        print*,"minimum acceptable epicentral distance?"
        read*,distmin
        print*,"increment along ray path?"
        read*,d_incr
        print *,"station distribution mode?"
        read*,mode
        print *,"events per station?"
        read*,ipick



c--------------------------------------------2d parameterization
        nxtot=xtot/x_incr
        nytot=ytot/y_incr
        print*,"then ",nxtot," horizontal intervals, and "
     &,nytot," vertical intervals"
        print*,nxtot*nytot,"pixels"
        if(nxtot*nytot.gt.nmax)stop "too many"

	iseed=-1
        rtmp = ran2(iseed)

c-------------------------------------main loop starts here
        open(9,file="sources.txt")
        open(19,file="receivers.txt")
        open(29,file="paths.txt")

        xcen=xtot/2.
        ycen=ytot/2.

        krow=0
c
c start of main loop until krow == ndata 
c
c

 1      continue
c-------------------------generate a source, a receiver, and the straight path between them
        krow=krow+1
        write(chrow,"(i6.6)")krow!test
        if(krow.gt.ndata)goto5

c
c     generate ipick events per station
c
        if((krow.eq.1).or.(mod(real(krow),real(ipick)).eq.0))then
c           print *,krow,ipick
c
c     make new station
c    
           if(mode.eq.1)then
c     uniformly distributed
              xrran=xtot*ran2(iseed)
              yrran=ytot*ran2(iseed)
           else if(mode.eq.2)then
c     gaussian distributed around center
              xrran = -1.
              do while((xrran.lt.0).or.(xrran.gt.xtot))
                 xrran=xcen + xtot/7*gasdev(iseed)
              end do
              yrran = -1. 
              do while((yrran.lt.0).or.(yrran.gt.ytot))
                 yrran=ycen + ytot/7*gasdev(iseed)
              end do
           else if((mode.eq.3).or.(mode.eq.4))then
c     distributed around sides
              dist = 0.
              xrran=-1.
              yrran=-1.
              do while((dist.lt.xtot/2.5).or.(xrran.lt.0).or.
     &             (xrran.gt.xtot).or.(yrran.lt.0).or.
     &             (yrran.gt.ytot))
                 xrran=xcen + xtot/7*gasdev(iseed)
                 yrran=ycen + ytot/7*gasdev(iseed)
                 dist=sqrt((xrran-xcen)**2+(yrran-ycen)**2)
              end do
           endif
           
        endif
c
c record receivers
c
        write(19,*)xrran,yrran !receivers
c
c
c     make new source, if needed, else use old
c
        
c
c     loop til we have a source at least distmin away
        epidist=0.
        ii=0
        do while((epidist.lt.distmin).and.(ii.lt.nsource))

           ii=ii+1
c     use an old source? only a few percent get reused
           if(ran2(iseed).lt.0.04)then
              xsran=source(ii,1)
              ysran=source(ii,2)
              epidist=sqrt((xsran-xrran)**2+(ysran-yrran)**2)
           endif
        enddo
        if((ii.ge.nsource).or.(epidist.lt.distmin))then
c     didn't find an old source
c
c     make new
c
           epidist=0.
           xsran = -1.
           ysran = -1.
           do while((epidist.lt.distmin).or.
     &          ((xsran.lt.0).or.
     &          (xsran.gt.xtot).or.(ysran.lt.0).or.
     &          (ysran.gt.ytot)))
c
c     pick a source
c     
              if(mode.eq.3)then ! ~Gaussian in center
                 xsran=xcen + xtot/3.*gasdev(iseed)
                 ysran=ycen + ytot/3.*gasdev(iseed)
              else if (mode.eq.4)then ! around center
                 dist = 0.
                 do while(dist.lt.xtot/2.5)
                    xsran=xcen + xtot/7*gasdev(iseed)
                    ysran=ycen + ytot/7*gasdev(iseed)
                    dist=sqrt((xsran-xcen)**2+(ysran-ycen)**2)
                 end do
              else
                 xsran=xtot*ran2(iseed)
                 ysran=ytot*ran2(iseed)
              endif
              epidist=sqrt((xsran-xrran)**2+(ysran-yrran)**2)
           enddo
           nsource=nsource+1
           source(nsource,1)=xsran
           source(nsource,2)=ysran
           write(9,*)xsran,ysran !sources
        endif

c
c     paths
c
        write(29,*)xsran,ysran
        write(29,*)xrran,yrran             !paths
        write(29,'(A)')">"


        azi=asin(abs(xsran-xrran)/epidist)
        dx=d_incr*sin(azi)
        dy=d_incr*cos(azi)

        if(xrran.lt.xsran)dx=-dx
        if(yrran.lt.ysran)dy=-dy

        d=0.
        k=1
        x(1)=xsran
        el(1)=ysran
 3      k=k+1
           x(k)=x(k-1)+dx
           el(k)=el(k-1)+dy
           d=d+d_incr
           if(d.gt.epidist)goto4
        goto3
 4      k=k-1
        x(k)=xrran
        el(k)=yrran
        np=k

c        open(39,file="ray."//chrow)
c        do k=1,np
c           write(39,*)k,x(k),el(k)
c        enddo
c        close(39)
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

c        print*,"sr=",krow,"ipx0=",ipx0

cTEST
c        open(44,file="inters."//chrow)

c------------loop over all points on the ray path
	do k=2,np
           nx=(x(k)/x_incr)+1
           ny=(el(k)/y_incr)+1

           if(inpx(nx,ny,nxtot,nytot).ne.ipx0)then ! crossed over to another pixel
              if(ny.ne.ny0)then !vertical intersection
                 yint=min(ny,ny0)*y_incr
                 if(ny.lt.ny0)yint=ny*y_incr!DO I NEED THIS LINE?
                 a=(x(k)-x0)/(el(k)-y0)
                 b=x0-a*y0
                 xint=a*yint+b !find x of intersection by linear interpolation
              elseif(nx.ne.nx0)then !horizontal intersection
                 xint=min(nx,nx0)*x_incr
                 a=(el(k)-y0)/(x(k)-x0)
                 b=y0-a*x0
                 yint=a*xint+b !find y of intersection by linear interpolation
              endif
              index=index+1
              ds(index)=sqrt(((xint-xint0)**2)+((yint-yint0)**2))
c              write(44,*)xint,yint,xint0,yint0,ipx0!test
              indnz(index)=ipx0
              xint0=xint
              yint0=yint
              nx0=nx
              ny0=ny
              ipx0=inpx(nx0,ny0,nxtot,nytot)
           endif
           x0=x(k)
           y0=el(k)
        enddo ! ray path is done
c----------------------------------pixel where receiver is
        if(x(np).ne.xint0)then
           index=index+1
           ds(index)=sqrt(((x(np)-xint0)**2)+((el(np)-yint0)**2))
           nx=(x(np)/x_incr)+1
           ny=(el(np)/y_incr)+1           
           indnz(index)=inpx(nx,ny,nxtot,nytot)
        endif
c	stop
c        close(44)!test
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
c        open(77,file="row."//chrow)
c        do iw=1,index
c           write(77,*)iw,indnz(iw),ds(iw)
c        enddo
c     close(77)

c        pause
c        stop


c-------------------------------------main loop ends here
	goto1
 5      continue

        close(9)
        close(19)
        close(29)

        end

           function inpx(nx,ny,nxtot,nytot)
           inpx=(ny-1)*nxtot+nx
           return
           end
