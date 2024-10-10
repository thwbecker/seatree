         subroutine f_tRes( t     , tdot  ,
     &                      x     , v     , idiagt,
     &                      trhs  , tlhsa , tlhsb )
c----------------------------------------------------------------------
c
c  S.D.K. 3/25/92 IMPLICIT version
c 
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
     &           x(ndof,*)    , v(ndof,*)    , idiagt(*)  ,
     &           trhs(*)      , tlhsa(*)     , tlhsb(*)   ,   
     &           eval(6)      ,
     &           shdx(4,5)    , shdy(4,5)    , det(5)
c
c
      common /temp1 / el_rhs(4) , blkmhu   , pg(4,5)  , 
     &                vl(8)     , uq(5)    , vq(5)    ,
     &                conv      , tmass    , adiff    , 
     &                unorm     , eta      , xse      ,
     &                blkdif    , uxse     , ueta     ,
     &                tlhsl(4,5), tl(4)    , xl(2,4)  ,
     &                m         , k        , igadd    
c
c
      do i = 1 , nEGnp
         trhs(i) = zero
      enddo
      do i = 1 , nsizet
	 tlhsa(i) = zero
         tlhsb(i) = zero
      enddo
c
c.... loop over the element blocks
c
      do 1000 ivel = 1 , numel
c
c..... localize the temperature ,tdot , velocity
c
      do n = 1,4
        vl(2*n-1) = v(1, ien(ivel,n) ) 
        vl(2*n)   = v(2, ien(ivel,n) ) 
        xl(1,n)   = x(1, ien(ivel,n) ) 
        xl(2,n)   = x(2, ien(ivel,n) ) 
#ifdef PICARD
        tl(n) = tbc( ien(ivel,n) ) 
#else
        tl(n) = t( ien(ivel,n) ) 
#endif
      enddo
c
c....form the the values of the temperature and velocity at intergration points
c
      do intp = 1,5
        uq(intp)    = vl(1) * shl(1,intp)  + vl(3) * shl(2,intp) 
     &              + vl(5) * shl(3,intp)  + vl(7) * shl(4,intp)
c    
        vq(intp)    = vl(2) * shl(1,intp)  + vl(4) * shl(2,intp) 
     &              + vl(6) * shl(3,intp)  + vl(8) * shl(4,intp)  
      enddo
c
c.... call the global shape function
c
        call genshg(shdx, shdy, det, xl, ivel)
c
c
c...for PG shape functions (see Brooks Thesis)
c
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
c
         uxse = abs(uq(5)*eval(1) + vq(5)*eval(2)) * eval(5)  
         ueta = abs(uq(5)*eval(3) + vq(5)*eval(4)) * eval(6)   
c
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
         unorm = uq(5) * uq(5) + vq(5) * vq(5)
         if( unorm .gt. 0.000001 ) then
           adiff = (uxse * xse + ueta * eta)/(two*unorm)
         else
           adiff = zero
         end if
c
c loop over integration points 
c
      do node = 1,4
        do intp = 1,5
c
c..... form pg shape function's for intergration point intp
c
          pg(node,intp) = shl(node,intp) + 
     &      adiff * (uq(intp)*shdx(node,intp)
     &            +  vq(intp)*shdy(node,intp))
c
          el_rhs(intp) = zero
c
        enddo
      enddo
c
c....  form Residue term at intergation point INTP (done in f_tlhs)
c   
      call f_tlhs( tlhsa,   tlhsb,     shl,
     &              shdx,    shdy,    
     &               lmt,    trhs,   idiagt, ivel,
     &               det)
c
c..... assemble this block's element residual
c
      do n = 1,4
        if(lmt(ivel,n) .ne. 0) 
     &  trhs(lmt(ivel,n)) = trhs(lmt(ivel,n)) + el_rhs(n)
      enddo
c      
c.... end loop over element blocks
c
1000  continue
c
      call unfact(tlhsa, tlhsb, idiagt, neqt)
      call unback(tlhsa, tlhsb, trhs, idiagt, neqt)
c
c.... adjust group assembled residual (trhs) for boundray conditions
c.... (i.e replace with boundray value)
c
      do n = 1 , nEGnp
          if(idt(n) .ne. 0 ) then
#ifdef PICARD
            t(n) = (t(n)+ trhs(idt(n)))/2.0
            tdot(n) = zero
#else
            tdot(n) = trhs(idt(n))
            t(n)    = t(n) + alpha * dt * tdot(n) 
#endif
          else
            t(n) = tbc(n)
            tdot(n) = zero
          end if
      enddo
c
c... wrap b.c.
c
      if (nwrap .gt. 0) then
c
        do n = 1, nwrap+1
          iwrap = nEGnp - nwrap-1 + n
          t(iwrap) = t(n)
          tdot(iwrap) = tdot(n)
        enddo
      end if
c
c....  return
c
      return
      end
