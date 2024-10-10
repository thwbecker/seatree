      subroutine stres( x , v, t )
c
      implicit double precision (a-h,o-z)
c
c.... remove above card for single-precision operation
c
c
#include<common.h>
#include<fe.h>
#include<stress.h>
#include<materials.h>
c
      dimension x(nsd,*) , v(ndof,*) , t(*) ,
     &         shdx(4,5) , shdy(4,5) , det(5)
c
c
      common /temp1 /  xl(2,4), vl(8)    , strtmp(6) ,  
     &                 volume , smass    , tmass(4)  , evisc(5)
c
      real*8 maxsecinv, secinv

      maxsecinv = -1.0d6
      do 50 n = 1, 6
        do 45 i = 1, numnp
          stress(n,i) = zero
          pmass(i) = zero
45      continue
50    continue

      do 1000 ivel = 1 , numel
c
c..... localize coordinates and velocity
c
        do n=1,nen
          xl(1, n) = x(1,ien(ivel,n) ) 
          xl(2, n) = x(2,ien(ivel,n) ) 
          vl(2*n-1) = v(1,ien(ivel,n) ) 
          vl(2*n  ) = v(2,ien(ivel,n) ) 
        enddo
c
c.... call the global shape function
c
        call genshg(shdx, shdy, det, xl, ivel)
c
c
c find velocity derivative at element centers
c
c scale by material parameters
c
        call rheol ( x , v , t , ivel , evisc)
c
        do n=1,6
          strtmp(n) = 0
        enddo
c
c for temperature depenedent problems, you need stresses at
c the integration points, not at the center, as is in some
c old versions of this routine.   SDK 7/15/08
c
        do intp=1,4
c  normal stress x
         strtmp(1) = strtmp(1) + two*evisc(intp)*
     &       (shdx(1,intp)*vl(1) + shdx(2,intp)*vl(3)
     &      + shdx(3,intp)*vl(5) + shdx(4,intp)*vl(7))*det(intp)
c  normal stress y
         strtmp(2) = strtmp(2) + two*evisc(intp)*
     &       (shdy(1,intp)*vl(2) + shdy(2,intp)*vl(4)
     &      + shdy(3,intp)*vl(6) + shdy(4,intp)*vl(8))*det(intp)
c  shear stress
         strtmp(3) = strtmp(3) + evisc(intp)*
     &  + 0.5*(shdx(1,intp)*vl(2) + shdx(2,intp)*vl(4)
     &     + shdx(3,intp)*vl(6) + shdx(4,intp)*vl(8)
     &     + shdy(1,intp)*vl(1) + shdy(2,intp)*vl(3)
     &     + shdy(3,intp)*vl(5) + shdy(4,intp)*vl(7))*det(intp)
c  pressure
        strtmp(4) = strtmp(4) +
     &       (shdx(1,intp)*vl(1) + shdx(2,intp)*vl(3)
     &      + shdx(3,intp)*vl(5) + shdx(4,intp)*vl(7))*det(intp)
     &      +(shdy(1,intp)*vl(2) + shdy(2,intp)*vl(4)
     &      + shdy(3,intp)*vl(6) + shdy(4,intp)*vl(8))*det(intp)
c
c end intp loop
c
        enddo
        volume=det(1)+det(2)+ det(3) + det(4)
	do n=1,6
	  strtmp(n) = strtmp(n)/volume
	enddo
        strtmp(4) =-evisc(5)*alam(mat(ivel))*strtmp(4) 
        strtmp(5) = strtmp(2) - strtmp(4) 
        strtmp(6) = evisc(5) 
c
c assemble element stress contributions to nodes
c
        do node = 1,4
c
c calculte element "weight"
c
        tmass(node) = det(1)*shl(node,1) + det(2)*shl(node,2)
     &              + det(3)*shl(node,3) + det(4)*shl(node,4)
        do n=1,6
          stress(n,ien(ivel,node)) = stress(n,ien(ivel,node)) 
     &                            +    tmass(node)*strtmp(n)
        enddo
        enddo
        smass=tmass(1) + tmass(2) + tmass(3) + tmass(4)
        volume=det(1)+det(2)+ det(3) + det(4)
c
        do node = 1,4
          tmass(node) = tmass(node) * volume / smass
          pmass(ien(ivel,node)) = pmass(ien(ivel,node))+tmass(node)
        enddo
c
c end loop over elements
c
1000  continue
c
c
      do n=1,6
        do i=1, numnp
          stress(n,i) = stress(n,i) / pmass(i)
        enddo
      enddo
c
      do n=1,6
        do i=1,nodebn
c Hughes for testing purposes
        stress(n,nb(1,i))=two*stress(n,nb(1,i))-stress(n,nb(2,i))
c234567890123456789012345678901234567890123456789012345678901234567890
        enddo
      enddo
c
c
c
      return
      end
      subroutine prtstr(x,stress)
c
c.... program to print stress data
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      include 'common.h'
c
      dimension x(ndof,*), stress(6,*)

c
      write(imout,1000) ndof,nelx,nelz,numnp,nstep,time
      write(imout,1500)
      do 100 n=1,numnp
      write(imout,2000) n,(x(i,n),i=1,2),(stress(j,n),j=1,6)
  100 continue
c
c
      return
c
 1000 format(5i10,f10.6)
 1500 format('  node   x1        x2           txx            tzz',
     &       '      txz         P    ')
 2000 format(1x,i7,1x,2(1pe13.7,1x),6(1pe12.5,1x))
      end
