         subroutine f_tRes( t     ,tdot   ,
     &                      x     , v     , idiag ,
     &                      tlhs  , trhs  )
c----------------------------------------------------------------------
c
c  This routine calculates the R.H.S. temperature Residual
c  and solves and predict's for the new temperature  
c----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
c
#include<common.h>
#include<materials.h>
#include<fe.h>
#include<bc.h>
c
c
      dimension  t(*)         , tdot(*)      ,
     &           x(ndof,*)    , v(ndof,*)    , idiag(*)   ,
     &           tlhs(*)      , trhs(*)      , eval(6)    ,
     &           shdx(4,5)    , shdy(4,5)    , det(5)
c
      real*8        el_rhs(4),blkmhu, pg(4,5) ,
     &                conv    ,tmass  , adiff ,
     &                unorm   ,eta   , xse   ,
     &                blkdif  ,uxse   , ueta , 
     &                vl(8)   , uq(5) , vq(5) ,
     &                txq(5)  , tyq(5), tq(5) ,
     &                tdq(5)  , tl(4) , tdl(4),
     &                xl(2,4)
c
      if (npass .gt. 1 ) then
        flag = one
      else 
        flag = zero
      end if
c
      do 100 i = 1 , nEGnp
         trhs(i) = zero
100   continue
c
c.... loop over the element blocks
c
      do 1000 ivel = 1 , numel
c..... localize the temperature ,tdot , velocity
c
        xl(1,1) = x(1, ien(ivel,1) ) 
        xl(2,1) = x(2, ien(ivel,1) ) 
        xl(1,2) = x(1, ien(ivel,2) ) 
        xl(2,2) = x(2, ien(ivel,2) ) 
        xl(1,3) = x(1, ien(ivel,3) ) 
        xl(2,3) = x(2, ien(ivel,3) ) 
        xl(1,4) = x(1, ien(ivel,4) ) 
        xl(2,4) = x(2, ien(ivel,4) ) 
        vl(1) = v(1, ien(ivel,1) ) 
        vl(2) = v(2, ien(ivel,1) ) 
        vl(3) = v(1, ien(ivel,2) ) 
        vl(4) = v(2, ien(ivel,2) ) 
        vl(5) = v(1, ien(ivel,3) ) 
        vl(6) = v(2, ien(ivel,3) ) 
        vl(7) = v(1, ien(ivel,4) ) 
        vl(8) = v(2, ien(ivel,4) ) 
        tl(1) = t( lmt(ivel,1) ) 
        tl(2) = t( lmt(ivel,2) ) 
        tl(3) = t( lmt(ivel,3) ) 
        tl(4) = t( lmt(ivel,4) )
        if( npass .gt. 1) then 
           tdl(1) = tdot(  lmt(ivel,1) ) 
           tdl(2) = tdot(  lmt(ivel,2) ) 
           tdl(3) = tdot(  lmt(ivel,3) ) 
           tdl(4) = tdot(  lmt(ivel,4) )
        endif
c
c.... call the global shape function
c
        call genshg (shdx, shdy , det , xl , ivel)
c
c
c....form the the values of the temperature and velocity at intergration points
c
        do intp = 1, 5
          uq(intp) = vl(1) * shl(1,intp)  + vl(3) * shl(2,intp) 
     &             + vl(5) * shl(3,intp)  + vl(7) * shl(4,intp)
c    
          vq(intp) = vl(2) * shl(1,intp)  + vl(4) * shl(2,intp) 
     &             + vl(6) * shl(3,intp)  + vl(8) * shl(4,intp)
c    
          txq(intp) = tl(1)*shdx(1,intp) + tl(2)*shdx(2,intp)
     &              + tl(3)*shdx(3,intp) + tl(4)*shdx(4,intp)
c
          tyq(intp) = tl(1)*shdy(1,intp) + tl(2)*shdy(2,intp)
     &              + tl(3)*shdy(3,intp) + tl(4)*shdy(4,intp)
c    
        enddo
c
        if ( npass .gt. 1 ) then
          do intp = 1, 4
            tdq(intp) = tdl(1) * shl(1,intp) + tdl(2) * shl(2,intp) 
     &                + tdl(3) * shl(3,intp) + tdl(4) * shl(4,intp)
          enddo           
        end if          
c
c.... calculate min element dimension per node
c
        exse1 = pt5*( xl(1,2) + xl(1,3) - xl(1,4) - xl(1,1) )
        exse2 = pt5*( xl(2,2) + xl(2,3) - xl(2,4) - xl(2,1) )
        eeta1 = pt5*( xl(1,3) + xl(1,4) - xl(1,1) - xl(1,2) )
        eeta2 = pt5*( xl(2,3) + xl(2,4) - xl(2,1) - xl(2,2) )
        hxse = sqrt(exse1*exse1 + exse2*exse2)
        heta = sqrt(eeta1*eeta1 + eeta2*eeta2)
        eval(1) = exse1/hxse
        eval(2) = exse2/hxse
        eval(3) = eeta1/heta
        eval(4) = eeta2/heta
        eval(5) = hxse
        eval(6) = heta
        uxse = abs(uq(5)*eval(1) + vq(5)*eval(2) ) * eval(5)  
        ueta = abs(uq(5)*eval(3) + vq(5)*eval(4) ) * eval(6)   
c
        blkmhu = dmhu(mat(ivel)) 
        blkdif = diff(mat(ivel))
c
        if( uxse .gt. two*blkdif)  then
          xse = one - two*blkdif/uxse
        else
          xse = zero
        end if
c
        if( ueta .gt. two*blkdif)  then
           eta = one - two*blkdif/ueta
        else
           eta = zero
        end if
c
c..... form pg shape function's for intergration point 1
c
        unorm = uq(5) * uq(5) + vq(5) * vq(5)
        if( unorm .gt. 0.000001 ) then
           adiff = (uxse * xse + ueta * eta)/(two*unorm)
        else
           adiff = zero
        end if
c       write(6,*) "adiff",ivel,adiff,unorm,uq(5),vq(5),
c    &     vl(1), vl(2), vl(3),vl(4)
c    &     pg(1,1), pg(2,2), pg(3,3), pg(4,4)
        do node = 1, 4
          do intp = 1, 4
            pg(node,intp) = shl(node,intp)  
     &        + adiff*(uq(intp)*shdx(node,intp)
     &        +        vq(intp)*shdy(node,intp))
          enddo
          el_rhs(node) = zero
        enddo
c 
c....  form the Residue term
c   
        do intp=1,4
          conv = uq(intp) * txq(intp) + vq(intp) * tyq(intp)
          tmass = tdq(intp) * flag + conv - blkmhu
          do node=1,4
            el_rhs(node) = el_rhs(node)
     &      -det(intp)*( tmass * pg(node,intp)
     &      +blkdif * (shdx(node,intp) * txq(intp) 
     &               + shdy(node,intp) * tyq(intp)))
          enddo
        enddo
c 
c..... assemble this block's element residual
c
c       write(6,*) "el_rhs ", ivel,
c    &     el_rhs(1), el_rhs(2), el_rhs(3), el_rhs(4)
        do node = 1,4
          trhs(lmt(ivel,node)) = trhs(lmt(ivel,node)) + el_rhs(node)
        enddo
c      
c.... end loop over element blocks
c
1000  continue
c
c.... adjust group assembled residual (trhs) for boundray conditions
c.... (i.e replace with boundray value)
c
      do n = 1 , nEGnp
        if(idt(n) .eq. 0 ) trhs(n) = zero
      enddo
c
c.... correct for new temperature
c
      if (npass .eq. 1) then
        do n = 1 , nEGnp
          tdot(n) = tlhs(n) * trhs(n)
          t(n)    = t(n) + alpha * dt * tdot(n) 
        enddo
      else
        do n = 1 , nEGnp
          tdot(n) = tdot(n) + tlhs(n) * trhs(n)
          t(n)    = t(n) + alpha * dt * tlhs(n) * trhs(n)
        enddo
      end if
c
c... adjust new temperature for boundary conditions
c
      do n = 1 , nEGnp
        if(idt(n) .eq. 0 ) t(n) = tbc(n)
c         write(6,*) n, t(n)
      enddo
c
      if (nwrap .gt. 0) then
        do n = 1, nwrap+1
          iwrap = nEGnp - nwrap-1 + n
          t(iwrap) = t(n)
          tdot(iwrap) = tdot(n)
        enddo
      endif
c
c....  return
c
      return
      end
