      subroutine  intemp( x ,  t , tdot )
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
#include<common.h>
#include<bc.h>
c
      dimension   x(nsd,1) , t(*) , tdot(*)
c
c
c
1     read (iin,*,err=1,end=999) pert , xsize , zsize 
      pix = 3.141592/xsize
      piz = 3.141592/zsize
      if ( inrstr .eq. 0 ) then
c
c.... use perturbed conductive start
c
        do 100 i = 1,numnp
          tdot(i) = zero
100     continue
        do 200 i = 1,numnp
          k = idt(i)
          if(k .ne. 0) then
      t(i) = (zsize-x(2,i))/zsize+pert*dcos(x(1,i)*pix)*dsin(x(2,i)*piz)
c     t(i) = 0.5d0+pert*dcos(x(1,i)*pix)*dsin(x(2,i)*piz)
          else
            t(i) = tbc(i)
          end if
200     continue
        time0=zero
        istep0=0
      else
c
c.... use restart file
c
        read(irsin,1000)istep0,time0
        read(irsin,*)
          do 301 n = 1, numnp
            read(irsin,1001) k,t(n),tdot(n)  
301       continue
      end if
1000  format(7x,i5,10x,f10.6)
1001  format(1x,i5,4(1pe15.8,1x))
c
c.... return
c
      return                           
999   call error ('init-tmp','end file',iin)
      end
