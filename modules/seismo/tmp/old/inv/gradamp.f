c-----defines the matrix corresponding to gradient damping
c-----and appends it to the A matrix, directly in memory.
c-----September 2006, I have found and fixed a bug in this routine--lapo
      subroutine gradamp4(weight,weightv,nnn,nelp,xxx,indx,mpoin,t,m,n,
     &     nonz,nlay,koffs,nsqrs,nsqtot,nlatzones,n_0,eq_incr,
     &     nlatzomax)
      dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
      dimension indx(nonz),xxx(nonz),mpoin(0:m),t(m)
      parameter(zer=0.)
      do 1000 l=1,nlay

         print*,l

         more=n_0*(L-1)+KOFFS
         vw=2.
         if((l.eq.1).or.(l.eq.nlay))vw=1.
         do i1=1,nlatzones
            rloin_2=180./nsqrs(i1)
            ifirst=nsqtot(i1)+1
            ilast=nsqtot(i1+1)
            do i2=ifirst,ilast
c--first, find ileft and iright
               ileft=i2-1
               if(i2.eq.ifirst)ileft=ilast
               iright=i2+1
               if(i2.eq.ilast)iright=ifirst
c--then, iup and idw...
               IF((1.LE.I2).AND.(I2.LE.nsqrs(1)))THEN
                  call coordsuper(nlatzomax,i2,rila,rilo,nsqrs,nlatzones,eq_incr)
                  dwla=rila-eq_incr
                  xlol=rilo-rloin_2
                  xlor=rilo+rloin_2
                  idwl=isqre(dwla,xlol,nsqrs,nsqtot,nlatzones,n_0,eq_incr)
                  idwr=isqre(dwla,xlor,nsqrs,nsqtot,nlatzones,n_0,eq_incr)
c--increment the matrix
                  xxx(nnn+1)=(3.+vw)*weight
                  indx(nnn+1)=i2+MORE
                  xxx(nnn+2)=-weight
                  indx(nnn+2)=iright+MORE
                  xxx(nnn+3)=-weight
                  indx(nnn+3)=ileft+MORE
                  nnn=nnn+3
c*****************bug fixed below, lapo 03.sep.2006
                  if(idwr.ge.idwl)then
                     do idw=idwl,idwr
                        xxx(nnn+1)=-weight/(idwr-idwl+1)
                        indx(nnn+1)=idw+MORE
                        nnn=nnn+1
                     enddo
                  elseif(idwr.lt.idwl)then
                     ndw=nsqtot(i1+2)-idwl+1
                     ndw=ndw+idwr-nsqtot(i1+1)
                     do idw=idwl,nsqtot(i1+2)
                        xxx(nnn+1)=-weight/ndw
                        indx(nnn+1)=idw+MORE
                        nnn=nnn+1
                     enddo
                     do idw=nsqtot(i1+1)+1,idwr
                        xxx(nnn+1)=-weight/ndw
                        indx(nnn+1)=idw+MORE
                        nnn=nnn+1
                     enddo		       
                  endif
c*****************
C====================================================
               elseIF((nsqrs(1)+1.LE.I2).AND.(I2.LE.n_0-nsqrs(nlatzones)))THEN
                  call coordsuper(nlatzomax,i2,rila,rilo,nsqrs,nlatzones,eq_incr)
                  dwla=rila-eq_incr
                  upla=rila+eq_incr
                  xlol=rilo-rloin_2
                  xlor=rilo+rloin_2
                  idwl=isqre(dwla,xlol,nsqrs,nsqtot,nlatzones,n_0,eq_incr)
                  idwr=isqre(dwla,xlor,nsqrs,nsqtot,nlatzones,n_0,eq_incr)
                  iupl=isqre(upla,xlol,nsqrs,nsqtot,nlatzones,n_0,eq_incr)
                  iupr=isqre(upla,xlor,nsqrs,nsqtot,nlatzones,n_0,eq_incr)
c--increment the matrix
                  xxx(nnn+1)=(4.+vw)*weight
                  indx(nnn+1)=i2+MORE
                  xxx(nnn+2)=-weight
                  indx(nnn+2)=iright+MORE
                  xxx(nnn+3)=-weight
                  indx(nnn+3)=ileft+MORE
                  nnn=nnn+3
c*****************bug fixed below, lapo 03.sep.2006
                  if(idwr.ge.idwl)then
                     do idw=idwl,idwr
                        xxx(nnn+1)=-weight/(idwr-idwl+1)
                        indx(nnn+1)=idw+MORE
                        nnn=nnn+1
                     enddo
                  elseif(idwr.lt.idwl)then
                     ndw=nsqtot(i1+2)-idwl+1
                     ndw=ndw+idwr-nsqtot(i1+1)
                     do idw=idwl,nsqtot(i1+2)
                        xxx(nnn+1)=-weight/ndw
                        indx(nnn+1)=idw+MORE
                        nnn=nnn+1
                     enddo
                     do idw=nsqtot(i1+1)+1,idwr
                        xxx(nnn+1)=-weight/ndw
                        indx(nnn+1)=idw+MORE
                        nnn=nnn+1
                     enddo		       
                  endif
c*****************
                  if(iupr.ge.iupl)then
                     do iup=iupl,iupr
                        xxx(nnn+1)=-weight/(iupr-iupl+1)
                        indx(nnn+1)=iup+MORE
                        nnn=nnn+1
                     enddo
                  elseif(iupr.lt.iupl)then
c     nup=nsqrs(i1-1)-iupl+1
c     nup=nup+iupl
                     nup=nsqtot(i1)-iupl+1
                     nup=nup+iupr-nsqtot(i1-1)
                     do iup=iupl,nsqtot(i1)
                        xxx(nnn+1)=-weight/nup
                        indx(nnn+1)=iup+MORE
                        nnn=nnn+1
                     enddo
                     do iup=nsqtot(i1-1)+1,iupr
                        xxx(nnn+1)=-weight/nup
                        indx(nnn+1)=iup+MORE
                        nnn=nnn+1
                     enddo
                  endif
c*****************
C======================================================
               elseIF((n_0-nsqrs(nlatzones)+1.LE.I2).AND.(I2.LE.n_0))THEN
                  call coordsuper(nlatzomax,i2,rila,rilo,nsqrs,nlatzones,eq_incr)
                  upla=rila+eq_incr
                  xlol=rilo-rloin_2
                  xlor=rilo+rloin_2
                  iupl=isqre(upla,xlol,nsqrs,nsqtot,nlatzones,n_0,eq_incr)
                  iupr=isqre(upla,xlor,nsqrs,nsqtot,nlatzones,n_0,eq_incr)
c--increment the matrix
                  xxx(nnn+1)=(3.+vw)*weight
                  indx(nnn+1)=i2+MORE
                  xxx(nnn+2)=-weight
                  indx(nnn+2)=iright+MORE
                  xxx(nnn+3)=-weight
                  indx(nnn+3)=ileft+MORE
                  nnn=nnn+3
c*****************
                  if(iupr.ge.iupl)then
                     do iup=iupl,iupr
                        xxx(nnn+1)=-weight/(iupr-iupl+1)
                        indx(nnn+1)=iup+MORE
                        nnn=nnn+1
                     enddo
                  elseif(iupr.lt.iupl)then
c     nup=nsqrs(i1-1)-iupl+1
c     nup=nup+iupl
                     nup=nsqtot(i1)-iupl+1
                     nup=nup+iupr-nsqtot(i1-1)
                     do iup=iupl,nsqtot(i1)
                        xxx(nnn+1)=-weight/nup
                        indx(nnn+1)=iup+MORE
                        nnn=nnn+1
                     enddo
                     do iup=nsqtot(i1-1)+1,iupr
                        xxx(nnn+1)=-weight/nup
                        indx(nnn+1)=iup+MORE
                        nnn=nnn+1
                     enddo
                  endif
c*****************
               ENDIF
c--vertical damping
               if(l.eq.1)then
                  xxx(nnn+1)=-weightV
                  indx(nnn+1)=i2+n_0+more
                  nnn=nnn+1
               elseif((l.gt.1).and.(l.lt.nlay))then
                  xxx(nnn+1)=-weightV
                  indx(nnn+1)=i2-n_0+more
                  xxx(nnn+2)=-weightV
                  indx(nnn+2)=i2+n_0+more
                  nnn=nnn+2
               elseif(l.eq.nlay)then
                  xxx(nnn+1)=-weightV
                  indx(nnn+1)=i2-n_0+more
                  nnn=nnn+1
               endif
               nelp=nelp+1
               mpoin(nelp)=nnn
               t(nelp)=zer
            enddo
         enddo
 1000 continue
      end
      
