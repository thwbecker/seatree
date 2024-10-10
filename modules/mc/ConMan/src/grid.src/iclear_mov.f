      subroutine iclear (iclr, n)
c
c----------------------------------------------------------------------
c
c  program to clear an integer array
c
c input:
c  n             : number of integers to be zeroed
c
c output:
c  iclr (n)      : the array to be zeroed
c----------------------------------------------------------------------
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension iclr(1)
c
      do 100 i = 1, n
        iclr(i) = 0
100   continue
c
      return
      end
      
      subroutine imove (idest, isourc, n)
c
c----------------------------------------------------------------------
c
c  program to move integer array
c
c input:
c  isourc (n)        : the array to be copied from
c  n                 : number of integers to be copied
c
c output:
c  idest  (n)        : the resulting array (to be copied to)
c
c----------------------------------------------------------------------
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension idest(1), isourc(1)
c
      do 100 i = 1, n
        idest(i) = isourc(i)
100   continue
c
      return
      end
