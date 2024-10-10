      subroutine vadd(vrhs , v)
c
      use bc

      implicit double precision (a-h,o-z)
c
      include 'common.h'
      dimension vrhs(*) , v(ndof,*) 
c
      do j = 1 , ndof
      do i = 1 , numnp
        if(idv(j,i) .ne. 0 ) then
           v(j,i) = vrhs(idv(j,i))
        else
           v(j,i) = vbc(j,i)
        end if
      enddo
      enddo
c
c   fix up v for wrap-around b.c. 
c 
      if (nwrap .ne. 0 ) then
        do  n=1,nwrap
          iwrap = numnp - nwrap+n
          v(1,iwrap) = v(1,n)
          v(2,iwrap) = v(2,n)
        enddo
      end if
      return
      end
