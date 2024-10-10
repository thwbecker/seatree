      subroutine smove (destin, source, n)
c
c----------------------------------------------------------------------
c
c  program to move floating point scalar to an array
c
c input:
c  source             : the scalar to be copied from
c  n                  : number of floating points to be copied
c
c output:
c  destin (n)         : the resulting array (to be copied to)
c
c----------------------------------------------------------------------
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
c
      dimension destin(1),source(1)
c
      do 100 i = 1, n
        destin(i) = source(i)
100   continue
c
      return
      end
