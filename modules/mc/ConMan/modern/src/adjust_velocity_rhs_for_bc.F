      subroutine adjust_velocity_rhs_for_bc
     &  (lmv, stiff, iel, vbcl, vbcr)

      implicit double precision (a-h,o-z)

      include 'common.h'

      dimension lmv(numel,8),stiff(8,8),vbcl(8),vbcr(*)

      do i = 1, 8     ! looping column
        if (lmv(iel,i) .ne. 0) then
          do j = 1, 8    ! looping row
            vbcr( lmv(iel,i) ) = vbcr( lmv(iel,i) )
     &                         - stiff(i,j) * vbcl(j)
          enddo
        endif
      enddo

c....  return

      return
      end
