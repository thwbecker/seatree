      subroutine dyna_P( x , v, t)

c calculation of total stress 
c Refer to 'Mechanics in the Earth and Environmental Sciences'
c (Middleton and Wilcock, page 299)
c and implementation of momentum equation of Changyeol's note

      use fe
      use mstres
      use materials

      implicit double precision (a-h,o-z)
c
      include 'common.h'
c
      dimension x(nsd,*) , v(ndof,*) , t(*), tlw(4)
c
      dimension xl(2,4), vl(2,4), strtmp(6) , rhsl(8)  , 
     &          tmass(4), evisc(5) , tl(4)  , rl(4) , tlr(4)

c    volume , smass

c.... initialization

      do 50 n = 1, 4
        do 45 i = 1, numnp
          stress(n,i) = zero
          pmass(i) = zero
45      continue
50    continue

      do 1000 iel = 1 , numel
c
c..... localize coordinates and velocity
c
        do n=1,nen
          xl(1,n) = x( 1,ien(iel,n) ) 
          xl(2,n) = x( 2,ien(iel,n) ) 
          vl(1,n) = v( 1,ien(iel,n) ) 
          vl(2,n) = v( 2,ien(iel,n) )
          rl(n)   = dens(ien(iel,n) )
          tl(n)   = t(   ien(iel,n) )
          tlr(n)  = tref(ien(iel,n) )
          tlw(n)  = tl(n) + tlr(n)
        enddo

        do i = 1, 6
           strtmp(i) = zero
        enddo

        call rheol ( xl , vl , tlw , iel , evisc)

        strtmp(4) = strtmp(4) + 
     &      ( shdx(iel,1,5)*vl(1,1) + shdx(iel,2,5)*vl(1,2)
     &      + shdx(iel,3,5)*vl(1,3) + shdx(iel,4,5)*vl(1,4) 
     &      + shdy(iel,1,5)*vl(2,1) + shdy(iel,2,5)*vl(2,2)
     &      + shdy(iel,3,5)*vl(2,3) + shdy(iel,4,5)*vl(2,4) )

c add compressible term to Vx,x, Vy,y and P' for the incompressible assumption
       
        bdens=(rl(1)*shl(1,5) + rl(2)*shl(2,5)
     &        +rl(3)*shl(3,5) + rl(4)*shl(4,5))
        uq = (vl(1,1)*shl(1,5) + vl(1,2)*shl(2,5)
     &     +  vl(1,3)*shl(3,5) + vl(1,4)*shl(4,5))
        vq = (vl(2,1)*shl(1,5) + vl(2,2)*shl(2,5)
     &     +  vl(2,3)*shl(3,5) + vl(2,4)*shl(4,5))
        rxq =(rl(1)*shdx(iel,1,5) + rl(2)*shdx(iel,2,5)
     &      + rl(3)*shdx(iel,3,5) + rl(4)*shdx(iel,4,5))
        ryq =(rl(1)*shdy(iel,1,5) + rl(2)*shdy(iel,2,5)
     &      + rl(3)*shdy(iel,3,5) + rl(4)*shdy(iel,4,5))

c terms for dynamic pressrue, P'
        strtmp(4) = strtmp(4) + (uq*rxq  + vq*ryq )/bdens 
c
c scale by material parameters
c
        strtmp(4) = evisc(5)*alam(mat(iel))*strtmp(4) 
c
c calculate element "weight"
c
        do n = 1, nen
          rhsl(n) = det(iel,1)*shl(n,1) + det(iel,2)*shl(n,2)
     &            + det(iel,3)*shl(n,3) + det(iel,4)*shl(n,4)
        enddo
c
c assemble element stress contribution to the node
c
        rhsl(1) = rhsl(1) * strtmp(4)
        rhsl(2) = rhsl(2) * strtmp(4)
        rhsl(3) = rhsl(3) * strtmp(4)
        rhsl(4) = rhsl(4) * strtmp(4)
        stress(4,ien(iel,1)) = stress(4,ien(iel,1)) + rhsl(1)
        stress(4,ien(iel,2)) = stress(4,ien(iel,2)) + rhsl(2)
        stress(4,ien(iel,3)) = stress(4,ien(iel,3)) + rhsl(3)
        stress(4,ien(iel,4)) = stress(4,ien(iel,4)) + rhsl(4)
c
c calculate pmass done above in fluxke left here for completeness
c
        do n = 1, nen
          tmass(n) = det(iel,1) * shl(n,1) + det(iel,2) * shl(n,2)
     &             + det(iel,3) * shl(n,3) + det(iel,4) * shl(n,4)
        enddo

        smass = tmass(1) + tmass(2) + tmass(3) + tmass(4)
        volume = det(iel,1) + det(iel,2) + det(iel,3) + det(iel,4)

        do n = 1, nen 
          tmass(n) = tmass(n) * volume / smass
        enddo
c
        do n = 1, nen
          pmass(ien(iel,n)) = pmass(ien(iel,n)) + tmass(n)
        enddo
c
c end loop over elements
c

1000  continue
c
c
c
        do i=1, numnp
          stress(4,i) = stress(4,i) / pmass(i)
        enddo
c
c stress treatment along the boudary and corner
c
      nbcst = 0

      if (nbcst .eq. 1) then
c
c Hughes for testing purposes
c
        do i=1,nodebn
          stress(4,nb(1,i))=two*stress(4,nb(1,i))-stress(4,nb(2,i))
        enddo

      elseif (nbcst .eq. 2) then
c
c Hughes for testing purposes (Noncorner, Boundary nodes, Hughes, 1987, 229p)
c

        do i=2, nelx
          stress(4,nb(1,i))=two*stress(4,nb(1,i))-stress(4,nb(2,i))   ! along bottom edge
        enddo

        do i=(nelx+3), (nodebn - 1) 
          stress(4,nb(1,i))=two*stress(4,nb(1,i))-stress(4,nb(2,i))   ! along top edge
        enddo

        do i= 2, nelz
          stress(4,i)=two*stress(4,i)-stress(4,(i+nelz+1))   ! along left edge
        enddo

        do i= (numnp-nelz+1), (numnp-1)
          stress(4,i)=two*stress(4,i)-stress(4,(i-nelz-1))   ! along right edge
        enddo
c
c  (External Corner nodes) 
c
       do ncorner = 1, 4
         if (ncorner .eq. 1 ) then     ! left-bottom corner
           nxa = 1
           nxb = 2
           nxc = nelz + 3
           nxd = nelz + 2
         endif  
         if (ncorner .eq. 2 ) then     ! right-bottom corner
           nxa = numnp - nelz
           nxb = numnp - nelz + 1
           nxc = numnp - 2 * nelz
           nxd = numnp - 2 * nelz - 1
         endif  
         if (ncorner .eq. 3 ) then     ! right-top corner
           nxa = numnp
           nxb = numnp - 1
           nxc = numnp - nelz - 2
           nxd = numnp - nelz - 1
         endif  
        if (ncorner .eq. 4 ) then     ! left-top corner
           nxa = nelz + 1
           nxb = nelz 
           nxc = 2 * nelz + 1
           nxd = 2 * (nelz + 1)
         endif  
c....  parameters
          x1a = x(1, nxa)     ;  x2a = x(2, nxa)
          x1b = x(1, nxb)     ;  x2b = x(2, nxb)
          x1c = x(1, nxc)     ;  x2c = x(2, nxc)
          x1d = x(1, nxd)     ;  x2d = x(2, nxd)
          pLB = x1c * x2d - x1d * x2c
          pLC = x1d * x2b - x1b * x2d
          pLD = x1b * x2c - x1c * x2b
          pLBbar = pLB + (x2c - x2d) * x1a + (x1d - x1c) * x2a
          pLCbar = pLC + (x2d - x2b) * x1a + (x1b - x1d) * x2a
          pLDbar = pLD + (x2b - x2c) * x1a + (x1c - x1b) * x2a  
          pL = pLB + pLC + pLD
              stress(4, nxa) = ( pLBbar * stress(4, nxb)
     &                               + pLCbar * stress(4, nxc)
     &                               + pLDbar * stress(4, nxd) ) / pL

       enddo

      endif ! end of stress treatment along boundary and corner

      return
      end
      
      
      
 
