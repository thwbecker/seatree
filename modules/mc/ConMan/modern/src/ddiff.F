      subroutine ddiff(x, t, diff)
c----------------------------------------------------------------------
c
c  This routine calculates equation of state for thermal expansivity variation
c  In the compressible case, expans will be multiplied at the front of 
c  Ra (Raliegh number).

c  input
c  x : node coodinate    (2,numnp) for 2D
c  t : node temperature (K) (numnp) for 2D
c  
c  output
c  diff: nodal thermal diffusivity   (numnp) for 2D
c----------------------------------------------------------------------
c
      use fe, only : ien! for ien
      use materials, only : tref
      implicit double precision (a-h,o-z)
c
      include 'common.h'
c
      dimension x(nsd,*), t(*), diff(*)
      dimension xl(2,4), tl(4), tlr(4), tlw(4)

c initialization
      do i = 1, 2
        do j = 1, 4
          xl(i,j) = zero ; tl(i) = zero
        enddo
      enddo
      
c.... loop over the element blocks
c
      do ivel = 1 , numel
c
c.... set up the parameters
c
        do n = 1, 4
          tl(n)       = t(  ien(ivel,n) )
          tlr(n)      = tref( ien(ivel,n) )
          tlw(n)      = tl(n) + tlr(n)
          xl(1,n)     = x(1, ien(ivel,n) )
          xl(2,n)     = x(2, ien(ivel,n) )
        enddo

c.... Equation of State for thermal conductivity
c.... thermal conductivity is non-dimensionalized value.  
c.... It is the divided value by a reference thermal conductivity at 298 K
c.... If the termal expansivity constant, use one as diff.

        do n = 1, 4
          diff(ien(ivel,n)) = one
        enddo
       
       enddo
c....     
      return
      end
      
      

