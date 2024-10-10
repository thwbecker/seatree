      subroutine stres( x , v, t )

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
      dimension x(nsd,*) , v(ndof,*) , t(*) ,
     &          tzzt(nelx+1) , tzzb(nelx+1) 
c
      dimension xl(2,4), vl(2,4), strtmp(6) , rhsl(8)  , 
     &          tmass(4), evisc(5) , tl(4)  , rl(4) , tlr(4), tlw(4)

c    volume , smass

c.... initialization
      do i = 1, 4
        rl(i) = zero ; tl(i) = zero ; tmass(i) = zero
      enddo

      do i = 1, 5
        evisc(i) = zero 
      enddo

      do i = 1, 6
        strtmp(i) = zero
      enddo

      do i = 1, nelx+1
        tzzt(i) = zero ; tzzb(i) = zero
      enddo

      do 50 n = 1, 5
        do 45 i = 1, numnp
          stress(n,i) = zero
          pmass(i) = zero
45      continue
50    continue

      if ( (itype .eq. 1) .or. (itype .eq. 2) ) then    ! ALA and TALA

      do iel = 1 , numel

c..... localize coordinates and velocity

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

        do i = 1, 6
          strtmp(i) = zero
        enddo
c
        call rheol ( xl , vl , tlw , iel , evisc)
c
c find velocity gradient at element centers 
c
        do intp = 1, nipt-1
c  normal gradient toward x (Vx,x)
          strtmp(1) = strtmp(1) + two * evisc(intp) * det(iel,intp) *
     &              (shdx(iel,1,intp)*vl(1,1)+shdx(iel,2,intp)*vl(1,2)
     &              +shdx(iel,3,intp)*vl(1,3)+shdx(iel,4,intp)*vl(1,4))
c  normal gradient toward y (Vy,y)
          strtmp(2) = strtmp(2) + two * evisc(intp) * det(iel,intp) * 
     &              (shdy(iel,1,intp)*vl(2,1)+shdy(iel,2,intp)*vl(2,2)
     &              +shdy(iel,3,intp)*vl(2,3)+shdy(iel,4,intp)*vl(2,4))
c  shear gradient toward x or y (Vx,y = Vy,x)
          strtmp(3) = strtmp(3) + evisc(intp) * det(iel,intp) * 
     &              (shdx(iel,1,intp)*vl(2,1)+shdx(iel,2,intp)*vl(2,2)
     &              +shdx(iel,3,intp)*vl(2,3)+shdx(iel,4,intp)*vl(2,4)
     &              +shdy(iel,1,intp)*vl(1,1)+shdy(iel,2,intp)*vl(1,2)
     &              +shdy(iel,3,intp)*vl(1,3)+shdy(iel,4,intp)*vl(1,4))
c sum of Vx,x and Vy,y for dynamic pressure term (P') based on the penalty method
          strtmp(4) = strtmp(4) + evisc(5) * det(iel,5) *
     &              (shdx(iel,1,5)*vl(1,1)+shdx(iel,2,5)*vl(1,2)
     &              +shdx(iel,3,5)*vl(1,3)+shdx(iel,4,5)*vl(1,4) 
     &              +shdy(iel,1,5)*vl(2,1)+shdy(iel,2,5)*vl(2,2)
     &              +shdy(iel,3,5)*vl(2,3)+shdy(iel,4,5)*vl(2,4) )
          strtmp(5) = strtmp(5) + evisc(intp) * det(iel,intp) 
c
c add compressible term to Vx,x, Vy,y and P' for the incompressible assumption
c 
          bdens= (rl(1)*shl(1,5) + rl(2)*shl(2,5)
     &         +  rl(3)*shl(3,5) + rl(4)*shl(4,5))*det(iel,5)
           vlx = (vl(1,1)*shl(1,5) + vl(1,2)*shl(2,5)
     &         +  vl(1,3)*shl(3,5) + vl(1,4)*shl(4,5))*det(iel,5)
           vly = (vl(2,1)*shl(1,5) + vl(2,2)*shl(2,5)
     &         +  vl(2,3)*shl(3,5) + vl(2,4)*shl(4,5))*det(iel,5)
           rxq = (rl(1)*shdx(iel,1,5)+rl(2)*shdx(iel,2,5)
     &         +  rl(3)*shdx(iel,3,5)+rl(4)*shdx(iel,4,5))*det(iel,5)
           ryq = (rl(1)*shdy(iel,1,5)+rl(2)*shdy(iel,2,5)
     &         +  rl(3)*shdy(iel,3,5)+rl(4)*shdy(iel,4,5))*det(iel,5)
c
c terms for dynamic pressrue, P'
c
           strtmp(4) = strtmp(4) 
     &               + evisc(5) * (vlx * rxq  + vly * ryq ) / bdens    
c
c terms caused by second coefficient of viscosity = - 2/3 * mhu
c
           strtmp(1)=strtmp(1) - two/three*evisc(intp)*det(iel,intp)*
     &              (shdx(iel,1,intp)*vl(1,1)+shdx(iel,2,intp)*vl(1,2)
     &              +shdx(iel,3,intp)*vl(1,3)+shdx(iel,4,intp)*vl(1,4)
     &              +shdy(iel,1,intp)*vl(2,1)+shdy(iel,2,intp)*vl(2,2)
     &              +shdy(iel,3,intp)*vl(2,3)+shdy(iel,4,intp)*vl(2,4))
           strtmp(2)=strtmp(2) - two/three*evisc(intp)*det(iel,intp)*
     &              (shdx(iel,1,intp)*vl(1,1)+shdx(iel,2,intp)*vl(1,2)
     &              +shdx(iel,3,intp)*vl(1,3)+shdx(iel,4,intp)*vl(1,4)
     &              +shdy(iel,1,intp)*vl(2,1)+shdy(iel,2,intp)*vl(2,2)
     &              +shdy(iel,3,intp)*vl(2,3)+shdy(iel,4,intp)*vl(2,4))
       enddo  ! end of intp

c scale by material parameters

        volume= det(iel,1) + det(iel,2) + det(iel,3) + det(iel,4)

        strtmp(1) = strtmp(1) / volume
        strtmp(2) = strtmp(2) / volume 
        strtmp(3) = strtmp(3) / volume 
        strtmp(4) = strtmp(4) * alam(mat(iel)) / volume 
        strtmp(5) = strtmp(5) / volume 
      
c calculate element "weight"

        do n = 1, 4
          rhsl(n) = det(iel,1) * shl(n,1) + det(iel,2) * shl(n,2)
     &            + det(iel,3) * shl(n,3) + det(iel,4) * shl(n,4)
        enddo

c assemble element stress contribution to the node
       
        do n=1,5
          rhsl(5) = rhsl(1) * strtmp(n)
          rhsl(6) = rhsl(2) * strtmp(n)
          rhsl(7) = rhsl(3) * strtmp(n)
          rhsl(8) = rhsl(4) * strtmp(n)
          stress(n,ien(iel,1)) = stress(n,ien(iel,1)) + rhsl(5)
          stress(n,ien(iel,2)) = stress(n,ien(iel,2)) + rhsl(6)
          stress(n,ien(iel,3)) = stress(n,ien(iel,3)) + rhsl(7)
          stress(n,ien(iel,4)) = stress(n,ien(iel,4)) + rhsl(8)
        enddo
c
c calculate pmass done above in fluxke left here for completeness
c
        do n = 1, 4
          tmass(n) = det(iel,1) * shl(n,1) + det(iel,2) * shl(n,2)
     &             + det(iel,3) * shl(n,3) + det(iel,4) * shl(n,4)
        enddo
        smass = tmass(1) + tmass(2) + tmass(3) + tmass(4)
        do n = 1, 4
          tmass(n) = tmass(n) * volume / smass
        enddo
        do n = 1, 4
          pmass(ien(iel,n)) = pmass(ien(iel,n)) + tmass(n)
        enddo
c
c end loop over elements
c
      enddo

      elseif ( (itype .eq. 3) .or. (itype .eq. 4) ) then  ! EBA and BA

      do iel = 1 , numel

c..... localize coordinates and velocity

        do n=1,nen
          xl(1,n) = x( 1,ien(iel,n) ) 
          xl(2,n) = x( 2,ien(iel,n) ) 
          vl(1,n) = v( 1,ien(iel,n) ) 
          vl(2,n) = v( 2,ien(iel,n) )
          tl(n)   = t(   ien(iel,n))
          tlr(n)  = tref(ien(iel,n))
          tlw(n)  = tl(n) + tlr(n)
        enddo

        do i = 1, 6
         strtmp(i) = zero
        enddo
c
        call rheol ( xl , vl , tlw , iel , evisc)
c
c find velocity gradient at element centers 
c
        do intp = 1, nipt-1
c
c  normal gradient toward x (Vx,x)
c
         strtmp(1) = strtmp(1) + two * evisc(intp) * det(iel,intp) *
     &             (shdx(iel,1,intp)*vl(1,1)+shdx(iel,2,intp)*vl(1,2)
     &             +shdx(iel,3,intp)*vl(1,3)+shdx(iel,4,intp)*vl(1,4))
c
c  normal gradient toward y (Vy,y)
c
         strtmp(2) = strtmp(2) + two * evisc(intp) * det(iel,intp) * 
     &             (shdy(iel,1,intp)*vl(2,1)+shdy(iel,2,intp)*vl(2,2)
     &             +shdy(iel,3,intp)*vl(2,3)+shdy(iel,4,intp)*vl(2,4))
c
c  shear gradient toward x or y (Vx,y = Vy,x)
c
         strtmp(3) = strtmp(3) + evisc(intp) * det(iel,intp) * 
     &             (shdx(iel,1,intp)*vl(2,1)+shdx(iel,2,intp)*vl(2,2)
     &             +shdx(iel,3,intp)*vl(2,3)+shdx(iel,4,intp)*vl(2,4)
     &             +shdy(iel,1,intp)*vl(1,1)+shdy(iel,2,intp)*vl(1,2)
     &             +shdy(iel,3,intp)*vl(1,3)+shdy(iel,4,intp)*vl(1,4))
c
c sum of Vx,x and Vy,y for dynamic pressure term (P') based on the penalty method
c
         strtmp(4) = strtmp(4) + evisc(5) * det(iel,5) *
     &             ( shdx(iel,1,5)*vl(1,1) + shdx(iel,2,5)*vl(1,2)
     &             + shdx(iel,3,5)*vl(1,3) + shdx(iel,4,5)*vl(1,4) 
     &             + shdy(iel,1,5)*vl(2,1) + shdy(iel,2,5)*vl(2,2)
     &             + shdy(iel,3,5)*vl(2,3) + shdy(iel,4,5)*vl(2,4) )

         strtmp(5) = strtmp(5) + evisc(intp) * det(iel,intp) 

        enddo  ! end of intp

c scale by material parameters

        volume= det(iel,1) + det(iel,2) + det(iel,3) + det(iel,4)
        strtmp(1) = strtmp(1) / volume
        strtmp(2) = strtmp(2) / volume 
        strtmp(3) = strtmp(3) / volume 
        strtmp(4) = strtmp(4) * alam(mat(iel)) / volume 
        strtmp(5) = strtmp(5) / volume 
c
c calculate element "weight"
c
        do n = 1, 4
          rhsl(n) = det(iel,1) * shl(n,1) + det(iel,2) * shl(n,2)
     &            + det(iel,3) * shl(n,3) + det(iel,4) * shl(n,4)
        enddo
c
c assemble element stress contribution to the node
c
        do n=1,5
          rhsl(5) = rhsl(1) * strtmp(n)
          rhsl(6) = rhsl(2) * strtmp(n)
          rhsl(7) = rhsl(3) * strtmp(n)
          rhsl(8) = rhsl(4) * strtmp(n)
          stress(n,ien(iel,1)) = stress(n,ien(iel,1)) + rhsl(5)
          stress(n,ien(iel,2)) = stress(n,ien(iel,2)) + rhsl(6)
          stress(n,ien(iel,3)) = stress(n,ien(iel,3)) + rhsl(7)
          stress(n,ien(iel,4)) = stress(n,ien(iel,4)) + rhsl(8)
        enddo
c
c calculate pmass done above in fluxke left here for completeness
c
        do n = 1, 4
          tmass(n) = det(iel,1) * shl(n,1) + det(iel,2) * shl(n,2)
     &             + det(iel,3) * shl(n,3) + det(iel,4) * shl(n,4)
        enddo
        smass = tmass(1) + tmass(2) + tmass(3) + tmass(4)
        volume= det(iel,1) + det(iel,2) + det(iel,3) + det(iel,4)
        do n = 1, 4
          tmass(n) = tmass(n) * volume / smass
        enddo
        do n = 1, 4
          pmass(ien(iel,n)) = pmass(ien(iel,n)) + tmass(n)
        enddo
c
c end loop over elements
c
      enddo

      endif
      do n=1,5
        do i=1, numnp
          stress(n,i) = stress(n,i) / pmass(i)
        enddo
      enddo

c stress treatment along the boudary and corner

      nbcst = 1

      if (nbcst .eq. 1) then

      do n=1,5
        do i=1,nodebn
          stress(n,nb(1,i))=two*stress(n,nb(1,i))-stress(n,nb(2,i))
        enddo
      enddo

      elseif (nbcst .eq. 2) then

c Hughes for testing purposes (Noncorner, Boundary nodes, Hughes, 1987, 229p)

       do n=1,5
        do i=2, nelx
          stress(n,nb(1,i))=two*stress(n,nb(1,i))-stress(n,nb(2,i))   ! along bottom edge
        enddo
       enddo

       do n=1,5
        do i=(nelx+3), (nodebn - 1) 
          stress(n,nb(1,i))=two*stress(n,nb(1,i))-stress(n,nb(2,i))   ! along top edge
        enddo
       enddo

       do n=1,5
        do i= 2, nelz
          stress(n,i)=two*stress(n,i)-stress(n,(i+nelz+1))   ! along left edge
        enddo
       enddo

       do n=1,5
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
          do n = 1, nipt
              stress(n, nxa) = ( pLBbar * stress(n, nxb)
     &                               + pLCbar * stress(n, nxc)
     &                               + pLDbar * stress(n, nxd) ) / pL

          enddo
        enddo

      endif ! end of stress treatment along boundary and corner
c
c calculation of geoid 
c      

c store Tau zz for the top surface
      do i = 1 ,(nelx+1) 
        mpoint  = nelz + 1 + (nelz + 1)*(i - 1)
        tzzt(i) = stress(2, mpoint) + stress(4,mpoint)
      enddo

c store Tau zz for the bottom surface
      do i = 1 ,(nelx+1)
        npoint  = 1 + (nelz + 1)*(i - 1)
        tzzb(i) = stress(2, npoint) + stress(4,npoint)
      enddo
      
      call geoid(x, t, tzzt, tzzb)
c
c
      return
      end
      
      
      
      subroutine prtstr(x)
c
c.... program to print stress data
c
      use mstres
      implicit double precision (a-h,o-z)
c
      include 'common.h'
c
      dimension x(ndof,*)

c
      write(imout,1000) ndof,nelx,nelz,numnp,time
      write(imout,1500)
      do 100 n=1,numnp
      write(imout,2000) n,(x(i,n),i=1,2),(stress(j,n),j=1,5)
  100 continue
c
c
      return
c
 1000 format(4i10,f10.6)
 1500 format('__node______x1__________x2____________txx             
     &tzz___________txz___________P_________Visc ')
 2000 format(1x,i7,1x,2(1pe13.7,1x),5(1pe12.5,1x))
      end


