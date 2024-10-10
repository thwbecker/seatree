      function coldot(a,b,n)
c
c.... program to compute the dot product of vectors stored column-wise
c
      implicit double precision (a-h,o-z)
c
      dimension a(*), b(*)
      common /const /zero  , pt25,    pt33,   pt5,     pt66,   one,
     &               onept5, two ,   three,  four,   sixten,  eps7
c
      coldot = zero
      psum = zero
c
!$omp parallel do if (n .ge. 20) reduction (+:psum)
      do i=1,n
      psum = psum + a(i)*b(i)
      enddo
!$omp end parallel do
      coldot = psum
c
      return
      end
