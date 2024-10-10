      subroutine adjbc(lmv, stiff, ivel, vbcl, vbcr)
c
c
      implicit double precision (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      include 'common.h'
c
      dimension lmv(numel,*),stiff(*),vbcl(*),vbcr(*)
      integer idiag(8)
c
      idiag(1) = 1
      idiag(2) = 3
      idiag(3) = 6
      idiag(4) = 10
      idiag(5) = 15
      idiag(6) = 21
      idiag(7) = 28
      idiag(8) = 36
c
      do j = 1, 8
      do i = 1, 8
        if ( lmv(ivel,i) .ne. 0 ) then
c
          if ( j .ge. i) then
            igadd = idiag( j ) - j + i
          else
            igadd = idiag( i ) - i + j
          endif
          vbcr( lmv(ivel,i) ) = vbcr( lmv(ivel,i) )
     &                        - stiff(igadd) * vbcl(j) 
c
        end if 
      enddo
      enddo
c
c....  return
c
      return
      end
