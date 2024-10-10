      subroutine tmass (x , tlhs )
c
c-------------------------------------------------------------------------
c
c   This routine calculates the contribution of the mass  for the thermal
c   left hand side
c
c-------------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
c.... remove above card for single-precision operation
c
c
#include<common.h>
#include<fe.h>
#include<bc.h>
c
        dimension  x(2,*), tlhs(*)  ,
     &          shdx(4,5), shdy(4,5), det(5)
c
    
c
        common /temp1 / tmassl(4) , volume , xl(2,4)
c
c.... loop over the elements blocks
c
      do ivel = 1,numel
        do n = 1,4
          xl(1,n) = x(1, ien(ivel,n) )
          xl(2,n) = x(2, ien(ivel,n) )
        enddo
c
c
c.... call the global shape function
c
        call genshg (shdx, shdy, det, xl, ivel)
c
         volume=det(1)+det(2)+det(3)+det(4)
c
        do n = 1, 4
          tmassl(n) =shl(n,1) * det(1) + shl(n,2) * det(2)
     &             + shl(n,3) * det(3) + shl(n,4) * det(4)
        enddo
c
        volume=volume/
     &           (tmassl(1)+tmassl(2)+tmassl(3)+tmassl(4))
c
        do n = 1, 4
          tlhs(ien(ivel,n)) = tlhs(ien(ivel,n)) + tmassl(n) * volume 
        enddo
c
c.... end element loop 
c
      enddo
c
c....  invert mass 
c
        nloop = nEGnp
        if (nwrap .gt. 0) nloop = nEGnp - nwrap - 1
        do 300 i = 1 , nloop
           k = idt(i)
           if(k.eq.0) then
              tlhs(i) = one
           else
              tlhs(i) = one/tlhs(i)
           end if
300   continue
c....  return
c
      return
      end
