        subroutine genshp (shldx,shldy)
c
c----------------------------------------------------------------------
c
c  This program generates the shape functions for bi-linear,and
c  calculates the min element dimension per node
c input:
c  x      (nsd,numnp)      : coordinates
c  ien    (numel,nen)      : ien array
c  numel                   : number of elements
c
c output:
c
c        shl  (nen,nipt)
c        shdx (numel,nen,nipt)
c        shdy (numel,nen,nipt)
c        det  (numel,nipt)   (X)
c        eval (numel,6)      (X)
c
c Note: the last five arrays are setup with element as the first index.
c       This should facilitate vectorization.
c
c----------------------------------------------------------------------
c
      use fe

      implicit double precision (a-h,o-z)
c
      include 'common.h'
c
      dimension sa(4), ta(4), sg(5), tg(5)  
      dimension shldx(nen,nipt), shldy(nen,nipt)

      do i = 1, nen
        sa(i) = zero ; ta(i) = zero
      enddo

      do i = 1, nipt
        sg(i) = zero ; tg(i) = zero
      enddo 
c
c.... set up parameters
c
      sa(1) = -pt5
      sa(2) =  pt5
      sa(3) =  pt5
      sa(4) = -pt5
c      
      ta(1) = -pt5
      ta(2) = -pt5
      ta(3) =  pt5
      ta(4) =  pt5 
c
      guass = one / dsqrt(three)
c
      sg(1) = -one * guass
      sg(2) =  one * guass
      sg(3) =  one * guass
      sg(4) = -one * guass
      sg(5) = zero * guass
c
      tg(1) = -one * guass
      tg(2) = -one * guass
      tg(3) =  one * guass
      tg(4) =  one * guass
      tg(5) = zero * guass 
c       
c   generate the generic shape functions
c
      do intp = 1,nipt
        do node = 1,nen
          shl(node,intp)   = (pt5 + sa(node) * sg(intp) ) 
     &                     * (pt5 + ta(node) * tg(intp) )
          shldx(node,intp) = sa(node)*(pt5 + ta(node) * tg(intp))
          shldy(node,intp) = ta(node)*(pt5 + sa(node) * sg(intp))
        enddo
      enddo
c
c.... return
c
      return
      end
