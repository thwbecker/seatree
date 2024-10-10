      subroutine form_velocity_rhs(x, t, vrhs, vbcr) 
c
      use fe
      use materials
      use mstres

      implicit double precision (a-h,o-z)
c
      include 'common.h'
c
      dimension x(nsd,*) , t(*) , vrhs(*) , vbcr(*), bdens(5)    
      dimension tq(5), tl(4), xl(2,4), rl(4), pq(5), plq(5)
c    
c
c---------------------------------------------------------------------------
c
c     This routinue calculates the R.H.S velocity Residual
c     
c---------------------------------------------------------------------------
c
c.... initialization
c
      do i = 1 , neqv
         vrhs(i) = vbcr(i)
      enddo


c.... loop over the element 

      do iel = 1 , numel

c..... localize the coordinate, temperature and velocity

        do n = 1 , nen
          xl(1,n) = x(1, ien(iel,n) )
          xl(2,n) = x(2, ien(iel,n) )
          rl(n)   = dens(ien(iel,n))
          tl(n)   = t( ien(iel,n) )
          pq(n)   = stress(4, ien(iel,n)) ! only needed for ALA case
        enddo

        do intp = 1, nipt  
          bdens(intp) = rl(1) * shl(1,intp) + rl(2) * shl(2,intp)   ! density
     &                + rl(3) * shl(3,intp) + rl(4) * shl(4,intp)
          tq(intp) = tl(1) * shl(1,intp) + tl(2) * shl(2,intp)      ! potential temperature
     &             + tl(3) * shl(3,intp) + tl(4) * shl(4,intp)
          plq(intp)= pq(1) * shl(1,intp) + pq(2) * shl(2,intp)      ! dynamic pressure
     &             + pq(3) * shl(3,intp) + pq(4) * shl(4,intp)
        enddo
c
c... form the element stiffness for this element
c
        blkra = ra(mat(iel))
        do node = 1 , nen
          if (lmv(iel,2*node  ) .ne. 0) then
            do intp = 1 , nipt-1
              vrhs(lmv(iel,2*node)) = vrhs(lmv(iel,2*node))
     &      + blkra*bdens(intp)*(tq(intp)*det(iel,intp)*shl(node,intp)) 
            enddo
          endif
        enddo
c
c pressure term
c
       if (itype .eq. 1) then    ! ALA
        do n = 1 , nen
          if (lmv(iel,2*n ) .ne. 0) then
c           do intp = 1 , 4 
                intp=5
                vrhs(lmv(iel,2*n)) = vrhs(lmv(iel,2*n))
     &   +      4.0*Di/cgamma*plq(intp)*det(iel,intp)*shl(n,intp)
c           enddo
          endif
        enddo
       endif 

c
c.... end loop over elements
c
      enddo
c
c....  return
c
      return
      end
