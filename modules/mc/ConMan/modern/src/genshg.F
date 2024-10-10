      subroutine  genshg ( x ) 
c----------------------------------------------------------------------
c
c  This program generates the globalshape functions for bi-linear,
c
c input:
c  xl     (nsd,nen)           : local coordinates
c  shldx  (nen ,nipt )        : local dx
c  shldy  (nen ,nipt )        : local dy
c  shl    (nen ,nipt )        : local shape functions
c
c output:
c
c        shdx (nen,nipt)
c        shdy (nen,nipt)
c        det  (nipt)
c
c Note: the last four arrays are setup with element as the first index.
c       This should facilitate vectorization.
c
c----------------------------------------------------------------------
c
c
      use fe

      implicit double precision (a-h,o-z)
c
      include 'common.h'
c       
      dimension xl(nsd, nen)
      dimension shldx(nen,nipt),shldy(nen,nipt),xs(2,2),
     &          x(2,numnp)
        
c
c.... loop over all the integration points
c
c.... initialization
      call genshp (shldx, shldy)

      do iel = 1, numel
        do n = 1, nen
         xl(1,n) = x(1,ien(numel,n))
         xl(2,n) = x(2,ien(numel,n))
        enddo
      do intp = 1 , nipt
c
c     find jocabian
c
         xs(1,1) = xl(1,1)*shldx(1,intp) + xl(1,2)*shldx(2,intp)
     &           + xl(1,3)*shldx(3,intp) + xl(1,4)*shldx(4,intp)
c
         xs(1,2) = xl(2,1)*shldx(1,intp) + xl(2,2)*shldx(2,intp)
     &           + xl(2,3)*shldx(3,intp) + xl(2,4)*shldx(4,intp)
c
         xs(2,1) = xl(1,1)*shldy(1,intp) + xl(1,2)*shldy(2,intp)
     &           + xl(1,3)*shldy(3,intp) + xl(1,4)*shldy(4,intp)
c
         xs(2,2) = xl(2,1)*shldy(1,intp) + xl(2,2)*shldy(2,intp)
     &           + xl(2,3)*shldy(3,intp) + xl(2,4)*shldy(4,intp)
c
c..... calculate the inverse jacobian
c
         det(iel,intp) = xs(1,1) * xs(2,2) - xs(1,2) * xs(2,1)
c
c.... check for zero determine
c
         if ( det(iel,intp) .le. zero )
     &      call error (' genshg  ','det-jacb',ivel)
c 
c.... continue inverse calculation
c.... and  find derivative with respect to global axes
c
        temp = dabs (one / det(iel,intp))
        do n = 1, nen
          shdx(iel,n,intp) = temp * ( xs(2,2) * shldx(n,intp) 
     &                          - xs(1,2) * shldy(n,intp) )
          shdy(iel,n,intp) = temp * (-xs(2,1) * shldx(n,intp) 
     &                          + xs(1,1) * shldy(n,intp) )
        enddo
c
c.... end of integration point loop
c
      enddo
      enddo
c
c.... return
c
      return
      end
