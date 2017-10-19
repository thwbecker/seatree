      include 'common_para.h'
      parameter(nmax=6720)      ! maximum number of free parameters

      dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
      dimension nsqrsh(nlatzhmax),nsqtoth(nlatzhmax+1)
      dimension iresolu(500000),inew(500000),iold(500000)
      dimension inewh(500000),ioldh(500000),ire2(500000)
      parameter(natamax=nmax*(nmax+1)/2)
      dimension gtgi(natamax),gtdi(nmax)
      dimension values(nmax)
      integer indx(nmax)
      character chn*6
	print*,"size of coarse pixels?"
	read*,eq_incr
	refgrid=eq_incr*1.
	nlatzones=180./eq_incr
	if(nlatzones.gt.nlatzomax)stop "coarse pixels too small"
	print*,"ratio to size of fine pixels?"
	read*,ifa
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
	if(nco+nfi.gt.nmax)then
	   print*,'but nmax=',n
	   pause
	endif
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

c      print*,"nfree=",nfree
        nfree=nco+nfi
      n=nfree
      if(nfree.gt.nmax)stop "not right"
      nz=0
      ndimata=nfree*(nfree+1)/2
c------------------------------reset
      do k=1,ndimata
         gtgi(k)=0.
      enddo
      do k=1,nfree
         gtdi(k)=0.
      enddo
      print*,"low-resolution parameters..."
      kpoin=0
      call rough_d_outside(gtgi,natamax,values,INDX,n,
     &     1.,kpoin,nz,nmax,gtdi,nsqrs,nsqtot
     &     ,nlatzones,eq_incr,n,iresolu,kireso,inew,
     &     EQ_HI,IFA,NSQRSH,NSQTOTH,NLATZOHI,INEWH,NUHI)
      print*,"high resolution parameters..."
      kpoin=0
      call rough_d_inside(gtgi,natamax,VALUES,INDX,n,1.,nz,nmax,
     &     nmax,gtdi,nsqrsh,nsqtoth,nlatzohi,eq_hi,nuhi,
     &     nsqrs,nsqtot,nlatzones,INEW,INEWH,ifa,eq_incr)
      print*,"store damping matrix on a file..."
      write(chn,'(i6.6)')nfree
      open(61,file='gtg.'//chn,access='direct',
     &     recl=4,form='unformatted')
      open(62,file='gtd.'//chn,access='direct',
     &     recl=4,form='unformatted')
      do kata=1,ndimata
         write(61,rec=kata)gtgi(kata)
      enddo
      do ib=1,nfree
         write(62,rec=ib)gtdi(ib)
      enddo
      close(61)
      close(62)
      print*,"...done"
      end

	subroutine rough_d_outside(ata,nata,g,indx,nblo,weight,nnn,
     &	nonz,m,atd,nsqrs,nsqtot,nlatzones,eq_incr,n,iresolu,kireso,inew,
     &	eq_hi,ifa,nsqrsh,nsqtoth,nlatzohi,inewh,numhi)
c--defines the matrix corresponding to gradient damping,
c--accounting for the surface gradient at each boundary between blocks.
	dimension g(nonz),atd(m),ata(nata)
	integer iresolu(500000),inew(500000),inewh(500000),indx(nonz)
	parameter(nlatzomax=180,nlatzhmax=720)
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
	dimension nsqrsh(nlatzhmax),nsqtoth(nlatzhmax+1)
cTEST
c	print*,numhi,nonz,m,eq_incr,nlatzones
c	pause

	kpoin0=kpoin
	nnorth=n/2
	indexch=(nlatzones/2)+1
	nsouth=nnorth+nsqrs(indexch)+1
	do i1=1,nlatzones
c--rloin_2=half the longitudinal size of the block
	   rloin_2=180./nsqrs(i1)
	   rloinhi_2=rloin_2/ifa
	   ifirst=nsqtot(i1)+1
	   ilast=nsqtot(i1+1)
	   do i2=ifirst,ilast
              nnn=0             !reset counter of nonzero row-entries
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
	            call coordsuper(i2,rila,rilo,nsqrs,nlatzones,eq_incr)
		    rloin=rilo+(rloin_2+rloinhi_2)
		    rlainstart=rila+(eq_incr-eq_hi)/2.
	            do kuks=1,ifa
		       rlain=rlainstart-((kuks-1)*eq_hi)
	               irihi=isqre(rlain,rloin,nsqrsh,nsqtoth,nlatzohi,numhi,eq_hi)
	               indx(nnn+1)=inewh(irihi)
	               g(nnn+1)=-weight/ifa
	               nnn=nnn+1
cTEST
c	print*,rlain,rloin
c	print*,nnn,iright,irihi,inewh(irihi),g(nnn)
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
c	      kpoin=kpoin+1
c	      mpoin(kpoin)=nnn
c	      rhs(kpoin)=0.
cTEST
c	write(*,*)(indx(kk),kk=1,nnn)
c	write(*,*)(g(kk),kk=1,nnn)
c        pause

              call contribution_ata(g,indx,ata,0.,atd,nonz,nata,nnn)
c----------------------------find block(s) to the south or north
              nnn=0             !reset counter of nonzero row-entries
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
c--define the row corresponding to variation wrt latitude (n hemisphere)
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
			  rladwhi=rila-(eq_incr+eq_hi)/2.
			  rlostart=xlomin+(fin_incr/2.)
	                  do kuks=1,ifa
	                     rlodwhi=rlostart+(kuks-1)*(fin_incr)
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

c	         kpoin=kpoin+1
c	         mpoin(kpoin)=nnn
c	         rhs(kpoin)=0.
cTEST
c	write(*,*)(indx(kk),kk=1,nnn)
c	write(*,*)(g(kk),kk=1,nnn)
c        pause

                 call contribution_ata(g,indx,ata,0.,atd,nonz,nata,nnn)

c==============================================
	      elseif((nsouth.le.i2).and.(i2.le.n))then
c--if i2 is in the s hemisphere (exclude the blocks bounded to the n
c--by the equator because the variation across the equator has already
c--been accounted for).
	         call coordsuper(i2,rila,rilo,nsqrs,nlatzones,eq_incr)
	         upla=rila+eq_incr
	         rilol=rilo-rloin_2
	         rilor=rilo+rloin_2
	         iupl=isqre(upla,rilol,nsqrs,nsqtot,nlatzones,n,eq_incr)
	         iupr=isqre(upla,rilor,nsqrs,nsqtot,nlatzones,n,eq_incr)
c--define the row corresponding to
c--variation wrt latitude (s hemisphere)
	         g(nnn+1)=weight
	         indx(nnn+1)=inew(i2)
	         nnn=nnn+1

	         peso=1./(iupr-iupl+1)
	         do iup=iupl,iupr
	            do k=1,kireso
	               if(iup.eq.iresolu(k))then
	                  call rang(idw,xlamin,xlamax,xlomin,xlomax,
     &	nsqrs,nsqtot,nlatzones,n,eq_incr)
c--block iup is in the high-resolution area; 
c--find the corresponding southernmost finer blocks:
	                  fin_incr=(360./nsqrs(i1-1))/ifa
	                  rlauphi=rila+((eq_incr+eq_hi)/2.)
			  rlostart=(xlomin+(fin_incr/2.))
	                  do kuks=1,ifa
	                     rlouphi=rlostart+(kuks-1)*(fin_incr)
	                     iuphi=isqre(rlauphi,rlouphi,nsqrsh,nsqtoth,nlatzohi,numhi,eq_hi)
	                     indx(nnn+1)=inewh(iuphi)
	                     g(nnn+1)=-weight*peso/ifa
	                     nnn=nnn+1
	                  enddo
	                  goto67
	               endif
	            enddo
	            g(nnn+1)=-weight*peso
	            indx(nnn+1)=inew(iup)
	            nnn=nnn+1
67	            continue
	         enddo
cTEST
c	write(*,*)(indx(kk),kk=1,nnn)
cTEST
c	write(*,*)(indx(kk),kk=1,nnn)
c	write(*,*)(g(kk),kk=1,nnn)
c       pause
                 call contribution_ata(g,indx,ata,0.,atd,nonz,nata,nnn)

c	         kpoin=kpoin+1
c	         mpoin(kpoin)=nnn
c	         rhs(kpoin)=0.
c==============================================
c==============================================
	      endif
47	      continue
	   enddo
	enddo
c	print*,'number of rows in the g.d.matrix:',kpoin-kpoin0
	return
	end
c--------------------------------------------------------------------
	subroutine rough_d_inside(ata,nata,g,indx,nblo,weight,nnn,nonz,
     &	m,atd,nsqrs,nsqtot,nlatzones,eq_incr,n,
     &	nsqrsco,nsqtotco,nlatzoco,inew,inewh,ifa,co_incr)

        dimension ata(nata)
	dimension indx(nonz),g(nonz),atd(m)
	parameter(nlatzomax=180,nlatzhmax=720)
	dimension nsqrsco(nlatzomax),nsqtotco(nlatzomax+1)
	dimension nsqrs(nlatzhmax),nsqtot(nlatzhmax+1)
	dimension inewh(500000),inew(500000)
cTEST
c	print*,n,co_incr,eq_incr,nlatzones,nlatzoco
c	pause

	kpoin0=kpoin
	print*,'gradamp=',weight

	nnorth=n/2
c	print*,'nnorth=',nnorth
	indexch=(nlatzones/2)+1
	nsouth=nnorth+nsqrs(indexch)+1
c	print*,'nsouth=',nsouth

c-----------------------------------------loop over latitudinal zones 
	do i1=1,nlatzones
	   rloin_2=180./nsqrs(i1)
	   ifirst=nsqtot(i1)+1
	   ilast=nsqtot(i1+1)
c--------------------------------loop over all blocks within lat zone
	   do i2=ifirst,ilast
              nnn=0             !reset counter of nonzero row entries
c----------------------E-W gradient
	    if(inewh(i2).eq.-1)goto64
	      iright=i2+1
	      if(i2.eq.ilast)iright=ifirst
	      if(inewh(iright).eq.-1)then
	       call coordsuper(iright,rila,rilo,nsqrs,nlatzones,eq_incr)
	       iright=isqre(rila,rilo,nsqrsco,nsqtotco,nlatzoco,n,co_incr)
	       indx(nnn+1)=inewh(i2)
	       indx(nnn+2)=inew(iright)
	       g(nnn+1)=weight/ifa
	       g(nnn+2)=-weight/ifa
	       nnn=nnn+2
	      else
	       indx(nnn+1)=inewh(i2)
	       indx(nnn+2)=inewh(iright)
	       g(nnn+1)=weight
	       g(nnn+2)=-weight
	       nnn=nnn+2
	      endif
cTEST
c	write(*,*)"EW",(indx(kk),g(kk),kk=1,nnn)
c	pause

              call contribution_ata(g,indx,ata,0.,atd,nonz,nata,nnn)
c	      kpoin=kpoin+1
c	      mpoin(kpoin)=nnn
c	      rhs(kpoin)=0.
c----------------------------------------------------------N-S gradient
c---------------------------------------------northern hemisphere
	      if((1.le.i2).and.(i2.le.nnorth))then
	         call coordsuper(i2,rila,rilo,nsqrs,nlatzones,eq_incr)
	         dwla=rila-eq_incr
		 rilol=rilo-rloin_2
		 rilor=rilo+rloin_2
	         idwla=dwla*100.
		 idwl=isqre(dwla,rilol,nsqrs,nsqtot,nlatzones,n,eq_incr)
		 idwr=isqre(dwla,rilor,nsqrs,nsqtot,nlatzones,n,eq_incr)
	         peso=1./(idwr-idwl+1)
	         iflig=0
c---------------------------------loop
	         do idw=idwl,idwr
                  nnn=0       !reset counter each time
	          if(inewh(idw).eq.-1)then
	           call coordsuper(idw,rila,rilo,nsqrs,nlatzones,eq_incr)
	           idw1=isqre(rila,rilo,nsqrsco,
     &	nsqtotco,nlatzoco,n,co_incr)
	           if(iflig.ne.idw1)then
	           g(nnn+1)=weight
	           indx(nnn+1)=inewh(i2)
	           nnn=nnn+1
	           indx(nnn+1)=inew(idw1)
	           g(nnn+1)=-weight
	           nnn=nnn+1
	           kpoin=kpoin+1
c	           mpoin(kpoin)=nnn
c	           rhs(kpoin)=0.

	           iflig=idw1
cTEST
c	write(*,*)"NS1",(indx(kk),kk=1,nnn)
                   call contribution_ata(g,indx,ata,0.,atd
     &	                ,nonz,nata,nnn)
	           endif
	          else
	           g(nnn+1)=weight*peso
	           indx(nnn+1)=inewh(i2)
	           nnn=nnn+1
	           g(nnn+1)=-weight*peso
	           indx(nnn+1)=inewh(idw)
	           nnn=nnn+1
c	           kpoin=kpoin+1
cTEST
c	write(*,*)"NS2",(indx(kk),kk=1,nnn)
                   call contribution_ata(g,indx,ata,0.,atd
     &	                ,nonz,nata,nnn)
c	           mpoin(kpoin)=nnn
c	           rhs(kpoin)=0.
	          endif
	         enddo
c------------------------------------------southern hemisphere
	      elseif((nsouth.le.i2).and.(i2.le.n))then
	         call coordsuper(i2,rila,rilo,nsqrs,nlatzones,eq_incr)
	         upla=rila+eq_incr
c	         ilol=(rilo-rloin_2)*100.+1
c	         ilor=(rilo+rloin_2)*100.-1
	         rilol=(rilo-rloin_2)
	         rilor=(rilo+rloin_2)
	         iupla=upla*100.
	         iupl=isqre(upla,rilol,nsqrs,nsqtot,nlatzones,n,eq_incr)
	         iupr=isqre(upla,rilor,nsqrs,nsqtot,nlatzones,n,eq_incr)
	         peso=1./(iupr-iupl+1)
	         iflig=0
	         do iup=iupl,iupr
                  nnn=0       !reset counter each time
	          if(inewh(iup).eq.-1)then
	           call coordsuper(iup,rila,rilo,nsqrs,nlatzones,eq_incr)
	           iup1=isqre(rila,rilo,nsqrsco,nsqtotco,nlatzoco,n,co_incr)
	           if(iflig.ne.iup1)then
	           g(nnn+1)=weight
	           indx(nnn+1)=inew(iup1)
	           nnn=nnn+1
	           indx(nnn+1)=inew(iup1)
	           g(nnn+1)=-weight
	           nnn=nnn+1
c	           kpoin=kpoin+1
c	           mpoin(kpoin)=nnn
c	           rhs(kpoin)=0.
cTEST
c	write(*,*)"NS3(S)",(indx(kk),kk=1,nnn)

                   call contribution_ata(g,indx,ata,0.,atd,nonz,nata,nnn)
	           iflig=iup1
	           endif
	          else
	           g(nnn+1)=weight*peso
	           indx(nnn+1)=inewh(i2)
	           nnn=nnn+1
	           g(nnn+1)=-weight*peso
	           indx(nnn+1)=inewh(iup)
	           nnn=nnn+1
c	           kpoin=kpoin+1
cTEST
c	write(*,*)"NS4(S)",(indx(kk),kk=1,nnn)
                   call contribution_ata(g,indx,ata,0.,atd,nonz,nata,nnn)
c	           mpoin(kpoin)=nnn
c	           rhs(kpoin)=0.
	          endif
	         enddo
	      endif
64	      continue
c--end of latitudinal and longitudinal loops:
	   enddo
cTEST
c	pause
	enddo
c	print*,'number of rows in the g.d.matrix:',kpoin-kpoin0
	return
	end
c--------------------------------------------------------------------
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
c**********************
	subroutine correspc(iresolu,kireso,n,ioffset,inew,iold,ico)
	dimension iresolu(500000),inew(500000),iold(500000)
	ico=0
	do i=1,n
	   do k=1,kireso
	   if(i.eq.iresolu(k))then
	      ico=ico+1
	      inew(i)=-1
	      goto42
	   endif
	   enddo
	   inew(i)=(i+ioffset)-ico
42	   iold(inew(i))=i
	enddo
	return
	end
c**********************
	subroutine corresph(iresolu,kireso,n,ioffset,inew,iold,ico)
	dimension iresolu(500000),inew(500000),iold(500000)
	ico=0
	do i=1,n
	   do k=1,kireso
	   if(i.eq.iresolu(k))then
	      inew(i)=(i+ioffset)-ico
	      goto42
	   endif
	   enddo
	   ico=ico+1
	   inew(i)=-1
42	   iold(inew(i))=i
	enddo
	return
	end
c**********************
       subroutine contribution_ata(row,index,ata,d,b,n,nata,nz)
c----given a row of A and corresponding datum augments AtA accordingly
        real*4 ata(nata),b(n),row(n),d
        integer index(n)
        do i=1,nz
          do j=i,nz
	    if(index(j).ge.index(i))then
               ind=(((index(j)-1)*index(j))/2)+index(i)
	    else
               ind=(((index(i)-1)*index(i))/2)+index(j)
	    endif
            ata(ind)=ata(ind)+row(i)*row(j)
          enddo
c--add to rhs vector the contribution of j-th row:
          b(index(i))=b(index(i))+(row(i)*d)
        enddo
        return
        end
c**********************
        subroutine contribution_a(row,d,n,io1,io2,io3,io4,
     &	irec,index,nz)
        real*4 row(n),d
        integer index(n)
        do i=1,nz
          irec=irec+1
          write(io1,rec=irec)row(i)
          write(io2,rec=irec)index(i)
        enddo
        write(io3,*)irec
        write(io4,*)d
        return
        end
c**********************
	subroutine coordsuper(nbloc,blocla,bloclo,nsqrs,nlatzones,eq_incr)
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
c**********************
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
