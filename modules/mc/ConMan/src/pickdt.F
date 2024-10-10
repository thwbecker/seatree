      subroutine pickdt( x , v )
c
c    routine to pick the courant time step
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
#include<common.h>
#include<materials.h>
#include<fe.h>
c
c
      dimension x(nsd,*) , v(ndof,*)
c
      common /temp1 /temp,  tempd, h,     blkdif,
     &               vl(8), uq(5), vq(5), xl(2,4)
c
c localize velocities
c
      uhmax=zero
      difmax=zero
      do ivel=1,numel
c
c.....form velocities at integration points
c
        do n = 1, 4
          vl(2*n-1) = v(1, ien(ivel,n) ) 
          vl(2*n  ) = v(2, ien(ivel,n) ) 
          xl(1,n  ) = x(1, ien(ivel,n) ) 
          xl(2,n  ) = x(2, ien(ivel,n) ) 
        enddo
c
        uq(5) = vl(1) * shl(1,5) + vl(3) * shl(2,5) 
     &        + vl(5) * shl(3,5) + vl(7) * shl(4,5)
c    
        vq(5) = vl(2) * shl(1,5) + vl(4) * shl(2,5) 
     &        + vl(6) * shl(3,5) + vl(8) * shl(4,5)
c  
c.... calculate min element dimension per node
c
        exse1 = pt5*( xl(1,2) + xl(1,3) - xl(1,4) - xl(1,1) )
        exse2 = pt5*( xl(2,2) + xl(2,3) - xl(2,4) - xl(2,1) )
        eeta1 = pt5*( xl(1,3) + xl(1,4) - xl(1,1) - xl(1,2) )
        eeta2 = pt5*( xl(2,3) + xl(2,4) - xl(2,1) - xl(2,2) )
        hxse = sqrt(exse1*exse1 + exse2*exse2)
        heta = sqrt(eeta1*eeta1 + eeta2*eeta2)
c
        blkdif=diff(mat(ivel))
        temp=dabs(uq(5))/hxse + dabs(vq(5))/heta
        h=dmin1(hxse,heta)
        tempd=blkdif/(h*h)

        if (temp .gt. uhmax) uhmax = temp
        if (tempd .gt. difmax) difmax = tempd
c
      enddo
c
      uhmax=dmax1(uhmax,difmax)
c
      dt=pt5/uhmax
c
      return
      end
