         subroutine form_temp_rhs_explicit
     &     (t, tdot, x, v, tlhs, trhs )
c----------------------------------------------------------------------
c
c  This routine calculates the R.H.S. temperature Residual
c  and solves and predict's for the new temperature  
c----------------------------------------------------------------------
c
      use materials
      use fe
      use bc

      implicit double precision (a-h,o-z)
c
      include 'common.h'
c
      dimension  t(*)         , tdot(*)      , 
     &           x(nsd,*)     , v(nsd,*)     , 
     &           tlhs(*)      , trhs(*)      , eval(6)    ,
     &           rl(4)        , diffl(4)     , tmassl(5)
c
      dimension el_rhs(4) , pg(4,5), uq(5)    , vq(5)     , 
     &                txq(5)   , tyq(5)      , tq(5)      ,
     &                tdq(5)   , tl(4)       , tdl(4)     ,
     &                xl(2,4)  , vl(2,4)     , uxq(5)     ,
     &                vyq(5)   , uyq(5)      , vxq(5)     ,
     &                tyqr(5)  , tlr(5)      , diffzq(5)  ,
     &                tlw(5)   , tqr(5)      , difflq(5)  ,
     &                rq(5)    , zq(5)       , evisc(5)
c
c.... initialization
c
      if (npass .gt. 1 ) then
        flag = one
      else 
        flag = zero
      end if
c
      do 100 i = 1 , numnp
         trhs(i) = zero
100   continue

c**************************************************************
c
c.... loop over the elements
c
      do 1000 iel = 1 , numel
c..... localize the temperature ,tdot , velocity
        do n=1, nen
        xl(1,n) = x(1, ien(iel,n) ) 
        xl(2,n) = x(2, ien(iel,n) ) 
        vl(1,n) = v(1, ien(iel,n) ) 
        vl(2,n) = v(2, ien(iel,n) )
        rl(n)   = dens(ien(iel,n) ) 
        diffl(n)= diff(ien(iel,n))
        tl(n)   = t( ien(iel,n) ) 
        tlr(n)  = tref( ien(iel,n) )

        tlw(n)  = tl(n) + tlr(n) 
        tdl(n) = zero
        if( npass .gt. 1) then 
           tdl(n) = tdot(  lmt(iel,n) ) 
        endif
      enddo

c....form the the values of the temperature and velocity at intergration points

        do intp = 1, 5
          uq(intp) = vl(1,1) * shl(1,intp) + vl(1,2) * shl(2,intp) 
     &             + vl(1,3) * shl(3,intp) + vl(1,4) * shl(4,intp)

          vq(intp) = vl(2,1) * shl(1,intp) + vl(2,2) * shl(2,intp) 
     &             + vl(2,3) * shl(3,intp) + vl(2,4) * shl(4,intp)

          zq(intp) = xl(2,1) * shl(1,intp) + xl(2,2) * shl(2,intp) 
     &             + xl(2,3) * shl(3,intp) + xl(2,4) * shl(4,intp)

          tq(intp) = tl(1) * shl(1,intp) + tl(2) * shl(2,intp)
     &             + tl(3) * shl(3,intp) + tl(4) * shl(4,intp)

          rq(intp) = rl(1) * shl(1,intp) + rl(2) * shl(2,intp)
     &             + rl(3) * shl(3,intp) + rl(4) * shl(4,intp)
   
          txq(intp)=tl(1)*shdx(iel,1,intp) + tl(2)*shdx(iel,2,intp)
     &             +tl(3)*shdx(iel,3,intp) + tl(4)*shdx(iel,4,intp)

          tyq(intp)=tl(1)*shdy(iel,1,intp) + tl(2)*shdy(iel,2,intp)
     &             +tl(3)*shdy(iel,3,intp) + tl(4)*shdy(iel,4,intp)

          uxq(intp)=vl(1,1)*shdx(iel,1,intp)+vl(1,2)*shdx(iel,2,intp)
     &             +vl(1,3)*shdx(iel,3,intp)+vl(1,4)*shdx(iel,4,intp)

          vyq(intp)=vl(2,1)*shdy(iel,1,intp)+vl(2,2)*shdy(iel,2,intp)
     &             +vl(2,3)*shdy(iel,3,intp)+vl(2,4)*shdy(iel,4,intp)

          uyq(intp)=vl(1,1)*shdy(iel,1,intp)+vl(1,2)*shdy(iel,2,intp)
     &             +vl(1,3)*shdy(iel,3,intp)+vl(1,4)*shdy(iel,4,intp)

          vxq(intp)=vl(2,1)*shdx(iel,1,intp)+vl(2,2)*shdx(iel,2,intp)
     &             +vl(2,3)*shdx(iel,3,intp)+vl(2,4)*shdx(iel,4,intp)
 
          tyqr(intp)=tlr(1)*shdy(iel,1,intp)+tlr(2)*shdy(iel,2,intp)
     &              +tlr(3)*shdy(iel,3,intp)+tlr(4)*shdy(iel,4,intp) 

          tqr(intp)= tlr(1) *shl(1,intp) + tlr(2) *shl(2,intp)
     &             + tlr(3) *shl(3,intp) + tlr(4) *shl(4,intp) 
     
          difflq(intp)=diffl(1)*shl(1,intp)+diffl(2)*shl(2,intp)
     &                +diffl(3)*shl(3,intp)+diffl(4)*shl(4,intp)
     
          diffzq(intp)=diffl(1)*shdy(iel,1,intp)
     &                +diffl(2)*shdy(iel,2,intp)
     &                +diffl(3)*shdy(iel,3,intp)
     &                +diffl(4)*shdy(iel,4,intp)
          
        enddo

        do intp = 1, nipt
          tdq(intp) = zero
        enddo           
        if ( npass .gt. 1 ) then
          do intp = 1, nipt-1
            tdq(intp) = tdl(1) * shl(1,intp) + tdl(2) * shl(2,intp) 
     &                + tdl(3) * shl(3,intp) + tdl(4) * shl(4,intp)
          enddo           
        end if          

c.... calculate min element dimension per node

        exse1 = pt5*( xl(1,2) + xl(1,3) - xl(1,4) - xl(1,1) )
        exse2 = pt5*( xl(2,2) + xl(2,3) - xl(2,4) - xl(2,1) )
        eeta1 = pt5*( xl(1,3) + xl(1,4) - xl(1,1) - xl(1,2) )
        eeta2 = pt5*( xl(2,3) + xl(2,4) - xl(2,1) - xl(2,2) )
        hxse = dsqrt(exse1*exse1 + exse2*exse2)
        heta = dsqrt(eeta1*eeta1 + eeta2*eeta2)
        eval(1) = exse1/hxse
        eval(2) = exse2/hxse
        eval(3) = eeta1/heta
        eval(4) = eeta2/heta
        eval(5) = hxse
        eval(6) = heta
        uxse = dabs(uq(5)*eval(1) + vq(5)*eval(2) ) * eval(5)  
        ueta = dabs(uq(5)*eval(3) + vq(5)*eval(4) ) * eval(6)   

        blkmhu = dmhu(mat(iel)) 
        blkdif = difflq(5)

        if( uxse .gt. two*blkdif)  then
          xse = one - two*blkdif/uxse
        else
          xse = zero
        end if

        if( ueta .gt. two*blkdif)  then
           eta = one - two*blkdif/ueta
        else
           eta = zero
        end if

c..... form pg shape function's for intergration point 1

        unorm = uq(5) * uq(5) + vq(5) * vq(5)
        if( unorm .gt. 0.000001 ) then
           adiff = (uxse * xse + ueta * eta)/(two*unorm)
        else
           adiff = zero
        end if

        do node = 1, nen
          do intp = 1, nipt-1
            pg(node,intp) = shl(node,intp)  
     &        + adiff*(uq(intp)*shdx(iel,node,intp)
     &        +        vq(intp)*shdy(iel,node,intp))
          enddo
          el_rhs(node) = zero
        enddo
 
c get evisc & det

        call rheol (xl, vl, tlw, iel, evisc)
c
        do intp = 1, nipt-1
c
c viscous dissipation term
c           
          Bvisc  = zero    
          if ((itype .eq. 3) .or. (itype .eq. 4)) then 
            Bvisc = two / three     ! EBA or BA
          endif
          phi = zero
          if (ra(mat(iel)) .gt. zero) then
            phi = (Di/ra(mat(iel))) * evisc(intp) *
     &      (two * (uxq(intp) * uxq(intp) + vyq(intp) * vyq(intp)) 
     &           + (vxq(intp) + uyq(intp))*(vxq(intp) + uyq(intp)) 
     &           + (Bvisc - two / three) 
     &           * (uxq(intp) + vyq(intp)) * (uxq(intp) + vyq(intp)) )
          endif
          if (itype .eq. 4) phi = zero
c
c adiabatic advection term
c
          adiab = Di*vq(intp)*tq(intp) 
          if (itype .eq. 3) then
            adiab = Di * vq(intp) * ( tq(intp) + T0/diff_T )
          elseif (itype .eq. 4) then
            adiab = zero
          endif
c
          tmassl(intp) =  rq(intp) * (tdq(intp) * flag - blkmhu
     &          + uq(intp)*txq(intp) + vq(intp)*tyq(intp)
     &          + adiab ) - phi 
        enddo
c end loop over intergration points   
c
        do node=1,nen
          do intp=1,nipt-1
            el_rhs(node) = el_rhs(node)
     &      -det(iel,intp)*( tmassl(intp) * pg(node,intp)
     &      + difflq(intp) * (shdx(iel,node,intp) * txq(intp) 
     &                      + shdy(iel,node,intp) * tyq(intp)))
          enddo
        enddo
        if ((itype .eq. 1) .or. (itype .eq. 2)) then
          do node=1,nen
            do intp=1,nipt-1
              el_rhs(node) = el_rhs(node)
     &          -det(iel,intp)*(
     &          diffzq(intp)*(tyq(intp)+tyqr(intp))+
     &          difflq(intp)*shdy(iel,node,intp)*tyqr(intp))
            enddo
          enddo
        endif
c 
c..... assemble this blocks element residual
c
        do node = 1,nen
          trhs(lmt(iel,node)) = trhs(lmt(iel,node)) + el_rhs(node)
        enddo
c      
c.... end loop over element blocks
c
1000  continue
ccccc

c.... adjust for heat flux boundary conditions.
c.... assumes user inputs heat-flux integrated across the boundary length 
c.... heat flux is applied on the element length unit
c.... For isuf(n,2), when n is 1:bottom, 2:right, 3:top and 4:left
c.... shl(n,5) = 0.25

      do n = 1, numsuf
        iel = isuf(n,1)
        if (isuf(n,2) .eq. 1) then
          trhs(lmt(iel,1))=trhs(lmt(iel,1))-suf(n,3)*two*shl(1,5)
          trhs(lmt(iel,2))=trhs(lmt(iel,2))-suf(n,3)*two*shl(2,5)
        end if
        if (isuf(n,2) .eq. 2) then
          trhs(lmt(iel,2))=trhs(lmt(iel,2))-suf(n,3)*two*shl(2,5)
          trhs(lmt(iel,3))=trhs(lmt(iel,3))-suf(n,3)*two*shl(3,5)
        end if
        if (isuf(n,2) .eq. 3) then
          trhs(lmt(iel,3))=trhs(lmt(iel,3))-suf(n,3)*two*shl(3,5)
          trhs(lmt(iel,4))=trhs(lmt(iel,4))-suf(n,3)*two*shl(4,5)
        end if
        if (isuf(n,2) .eq. 4) then
          trhs(lmt(iel,4))=trhs(lmt(iel,4))-suf(n,3)*two*shl(4,5)
          trhs(lmt(iel,1))=trhs(lmt(iel,1))-suf(n,3)*two*shl(1,5)
        end if
      enddo

c
c.... adjust group assembled residual (trhs) for boundray conditions
c.... (i.e replace with boundray value)
c
      do n = 1 , numnp
        if(idt(n) .eq. 0 ) trhs(n) = zero
      enddo
c
c.... correct for new temperature
c
      if (npass .eq. 1) then
        do n = 1 , numnp
          tdot(n) = tlhs(n) * trhs(n)
          t(n)    = t(n) + alpha * dt * tdot(n) 
        enddo
      else
        do n = 1 , numnp
          tdot(n) = tdot(n) + tlhs(n) * trhs(n)
          t(n)    = t(n) + alpha * dt * tlhs(n) * trhs(n)
        enddo
      end if
c
c... adjust new temperature for boundary conditions
c
      do n = 1 , numnp
        if(idt(n) .eq. 0 ) t(n) = tbc(n)
      enddo
c
      if (nwrap .gt. 0) then
cc        do n = 1, nwrap+1
        do n = 1, nwrap
cc          iwrap = numnp - nwrap-1 + n
          iwrap = numnp - nwrap + n
          t(iwrap) = t(n)
          tdot(iwrap) = tdot(n)
        enddo
      endif
c
c....  return
c
      return
      end


