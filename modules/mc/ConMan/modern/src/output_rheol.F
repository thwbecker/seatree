      subroutine output_rheol(x, v, t)

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
      dimension x(nsd,*) , v(ndof,*) , t(*) 
c
      dimension xl(2,4), vl(2,4), strtmp(6) , rhsl(8)  , 
     &          tmass(4), evisc(5) , tl(4)  , rl(4) , tlr(4), tlw(4)

c.... initialization
      do i = 1, 6
        strtmp(i) = zero
      enddo
      do i = 1, numnp
        stress(1,i) = zero
        pmass(i) = zero
      enddo

      do 1000 iel = 1 , numel
c..... localize coordinates and velocity

        write(6,*) "in output_rheol", iel
        do n=1,nen
          xl(1,n) = x( 1,ien(iel,n) ) 
          xl(2,n) = x( 2,ien(iel,n) ) 
          vl(1,n) = v( 1,ien(iel,n) ) 
          vl(2,n) = v( 2,ien(iel,n) )
          rl(n)   = dens(ien(iel,n))
          tl(n)   = t(   ien(iel,n))
          tlr(n)  = tref(ien(iel,n))
          tlw(n)  = tl(n) + tlr(n)
        enddo

        call rheol ( xl , vl , tl , iel , evisc)

        volume= det(iel,1) + det(iel,2) + det(iel,3) + det(iel,4)
       
        strtmp(1) = 0.0
        do intp = 1,4
          strtmp(1)=strtmp(1)+evisc(intp)*det(iel,intp)/volume
        enddo

c calculate element "weight"

        do n = 1, 4
          rhsl(n) = det(iel,1) * shl(n,1) + det(iel,2) * shl(n,2)
     &            + det(iel,3) * shl(n,3) + det(iel,4) * shl(n,4)
        enddo

c assemble element stress contribution to the node

        rhsl(5) = rhsl(1) * strtmp(1)
        rhsl(6) = rhsl(2) * strtmp(1)
        rhsl(7) = rhsl(3) * strtmp(1)
        rhsl(8) = rhsl(4) * strtmp(1)
        stress(1,ien(iel,1)) = stress(1,ien(iel,1)) + rhsl(5)
        stress(1,ien(iel,2)) = stress(1,ien(iel,2)) + rhsl(6)
        stress(1,ien(iel,3)) = stress(1,ien(iel,3)) + rhsl(7)
        stress(1,ien(iel,4)) = stress(1,ien(iel,4)) + rhsl(8)

c calculate pmass done above in fluxke left here for completeness

        do n = 1, 4
          tmass(n) = det(iel,1) * shl(n,1) + det(iel,2) * shl(n,2)
     &             + det(iel,3) * shl(n,3) + det(iel,4) * shl(n,4)
        enddo

        smass = tmass(1) + tmass(2) + tmass(3) + tmass(4)

c        volume= det(iel,1)   + det(iel,2)   + det(iel,3)   + det(iel,4)

        do n = 1, 4
          tmass(n) = tmass(n) * volume / smass
        enddo

        do n = 1, 4
          pmass(ien(iel,n)) = pmass(ien(iel,n)) + tmass(n)
        enddo

c end loop over elements
        
1000  continue

        do i=1, numnp
          stress(1,i) = stress(1,i) / pmass(i)
        enddo

c stress treatment along the boudary and corner

      nbcst = 1

      if (nbcst .eq. 1) then

      do n=1,4
        do i=1,nodebn
c Hughes for testing purposes
        stress(n,nb(1,i))=two*stress(n,nb(1,i))-stress(n,nb(2,i))
c234567890123456789012345678901234567890123456789012345678901234567890
        enddo
      enddo

      elseif (nbcst .eq. 2) then

c Hughes for testing purposes (Noncorner, Boundary nodes, Hughes, 1987, 229p)

       do n=1,4
        do i=2, nelx
          stress(n,nb(1,i))=two*stress(n,nb(1,i))-stress(n,nb(2,i))   ! along bottom edge
        enddo
       enddo

       do n=1,4
        do i=(nelx+3), (nodebn - 1) 
          stress(n,nb(1,i))=two*stress(n,nb(1,i))-stress(n,nb(2,i))   ! along top edge
        enddo
       enddo

       do n=1,4
        do i= 2, nelz
          stress(n,i)=two*stress(n,i)-stress(n,(i+nelz+1))   ! along left edge
        enddo
       enddo

       do n=1,4
        do i= (numnp-nelz+1), (numnp-1)
          stress(n,i)=two*stress(n,i)-stress(n,(i-nelz-1))   ! along right edge
        enddo
       enddo


c  (External Corner nodes) 
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
          do n = 1, 4
              stress(n, nxa) = ( pLBbar * stress(n, nxb)
     &                               + pLCbar * stress(n, nxc)
     &                               + pLDbar * stress(n, nxd) ) / pL

         enddo
       enddo

      endif ! end of stress treatment along boundary and corner

c     do i=1,numnp
c       write(42,*) x(1,i), x(2,i), stress(1,i)
c     enddo
c
      return
      end
