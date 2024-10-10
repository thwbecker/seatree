      subroutine f_vRes(x , t , vrhs , vbcr  ) 
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
#include<common.h>
#include<materials.h>
#include<fe.h>
c
      dimension x(2,*) , t(*) , vrhs(*) , vbcr(*),
     &           shdx(4,5)    , shdy(4,5)    , det(5)
c
c
c---------------------------------------------------------------------------
c
c     This routinue calculates the R.H.S velocity Residual
c     
c---------------------------------------------------------------------------
c
      dimension tq(5), tl(4), xl(2,4)
c
      do i = 1 , neqv
         vrhs(i) = vbcr(i)
      enddo
c
c.... loop over the element blocks
c
      do ivel = 1 , numel
c
c..... localize the coordinate, temperature
c
        do n = 1 , 4
          tl(n) = t( ien(ivel,n) ) 
          xl(1,n) = x(1, ien(ivel,n) )
          xl(2,n) = x(2, ien(ivel,n) )
        enddo
c
c....form the temperature at the integration points 
c
        do m = 1 , 5
          tq(m) = shl(1,m) * tl(1) + shl(2,m) * tl(2)
     &          + shl(3,m) * tl(3) + shl(4,m) * tl(4)
        enddo
c
c.... call the global shape function
c
        call genshg(shdx, shdy, det, xl, ivel)
c
c
c... form the element stiffness for this block
c
        do n = 1 , 4
          if (lmv(ivel,2*n) .ne. 0) then
            blkra = ra(mat(ivel))
            do m = 1 , 4
              vrhs(lmv(ivel,2*n)) = vrhs(lmv(ivel,2*n)) 
     &          + blkra * ( tq(m) * det(m) * shl(n,m) )
            enddo
          endif
        enddo
c
c.... end loop over element blocks
c
      enddo
c
c....  return
c
      return
      end
