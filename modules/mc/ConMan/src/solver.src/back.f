      subroutine back(a,b,idiag,neq)
c
c.... program to perform forward reduction and back substitution
c
      implicit double precision (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension a(*),b(*),idiag(*)
      common /const /zero  , pt25,    pt33,   pt5,     pt66,   one,    
     &               onept5, two ,   three,  four,   sixten,  eps7
c
c.... forward reduction
c
      jj = 0
c
      do 100 j=1,neq
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
      if (jcolht.gt.1)
     &   b(j) = b(j) - coldot(a(jjlast+1),b(j-jcolht+1),jcolht-1)
  100 continue
c
c.... diagonal scaling
c
      do 200 j=1,neq
      ajj = a(idiag(j))
      if (ajj.ne.zero) b(j) = b(j)/ajj
  200 continue
c
c.... back substitution
c
      if (neq.eq.1) return
      jjnext = idiag(neq)
c
      do 400 j=neq,2,-1
      jj     = jjnext
      jjnext = idiag(j-1)
      jcolht = jj - jjnext
      if (jcolht.gt.1) then
         bj = b(j)
         istart = j - jcolht + 1
         jtemp  = jjnext - istart + 1
c
         do 300 i=istart,j-1
         b(i) = b(i) - a(jtemp+i)*bj
  300    continue
c
      endif
c
  400 continue
c
      return                                                                    
      end                                                                       
