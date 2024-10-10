      subroutine rdflux( )
c
c reads surface force and flux cards
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
#include<common.h>
#include<bc.h>
c
      do 100 k=1,numsuf
1       read(iin,*,err=1,end=999) isuf(k,1), isuf(k,2), suf(k,1),
     &                            suf(k,2) , suf(k,3)
 100  continue
c
c.... return
c
      return
c
c.... end of file error handling
c	
 999  call error('rdflux  ','end file',iin)
c
      end
