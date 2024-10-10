      subroutine  genshg (shdx , shdy , det , xl , ivel ) 
c
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
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
#include<common.h>
#include<fe.h>
c       
      dimension xl(nsd,*)
      dimension shdx(nen,*) ,shdy(nen,*) ,det(*), xs(2,2)
c
c.... loop over all the integration points
c
      do intp = 1 , nipt
c
c     find jocabian

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
         det(intp) = xs(1,1) * xs(2,2) - xs(1,2) * xs(2,1)
c
c.... check for zero determine
c
         if ( det(intp) .le. zero )
     &      call error (' genshg  ','det-jacb',ivel)
c 
c.... continue inverse calculation
c.... and  find derivative with respect to global axes
c
        temp = one/det(intp)
        do n = 1, nen
          shdx(n,intp) = temp * ( xs(2,2) * shldx(n,intp) 
     &                          - xs(1,2) * shldy(n,intp) )
          shdy(n,intp) = temp * (-xs(2,1) * shldx(n,intp) 
     &                          + xs(1,1) * shldy(n,intp) )
        enddo
c
c.... end of integration point loop
c
      enddo
c
c.... return
c
      return
      end
