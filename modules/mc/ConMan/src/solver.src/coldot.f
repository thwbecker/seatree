      function coldot(a,b,n)
c
c.... program to compute the dot product of vectors stored column-wise
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
c
c
c.... remove above card for single-precision operation
c
      dimension a(*),b(*)
      common /const /zero  , pt25,    pt33,   pt5,     pt66,   one,
     &               onept5, two ,   three,  four,   sixten,  eps7
c
      zero = 0.0
      coldot = zero
c
      do 100 i=1,n
      coldot = coldot + a(i)*b(i)
  100 continue
c
      return
      end
