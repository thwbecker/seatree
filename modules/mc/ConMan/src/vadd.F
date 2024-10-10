       subroutine vadd(vrhs , v)
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
#include<common.h>
#include<bc.h>
       dimension vrhs(*) , v(*) 
c
c
       do 100 j = 1 , ndof
       do 100 i = 1 , numnp
          if(idv(j,i) .ne. 0 ) then
             v(j+ndof*(i-1)) = vrhs(idv(j,i))
          else
             v(j+ndof*(i-1)) = vbc(j,i)
          end if
100   continue
c
c   fix up v for wrap-around b.c. 
c 
      if (nwrap .ne. 0 ) then
        do 200 n=1,2*nwrap+2
          iwrap = nEGdf - 2*nwrap-2+n
          v(iwrap) = v(n)
200     continue
      end if
      return
      end
