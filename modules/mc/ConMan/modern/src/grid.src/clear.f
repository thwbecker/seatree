      subroutine clear (clr, n)
c
      implicit double precision (a-h,o-z)
c
c----------------------------------------------------------------------
c  program to clear a floating point array
c
c input:
c  n        : number of floating points to be zeroed
c
c output:
c  clr (n)  : the array to be zeroed
c----------------------------------------------------------------------
c
      include 'common.h'
c
      dimension clr(*)
c
      do 100 i = 1, n
        clr(i) = zero
100   continue
c
      return
      end
