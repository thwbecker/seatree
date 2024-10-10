      subroutine masflx()
c
c--------------------------------------------------------------------
c  subroutine for program conman for mean field and average nusselt
c  number calculation.
c--------------------------------------------------------------------
c
c  Added irow to argument list on 6/1/94 jjita
c
c
      use variables
      use fe
      use materials
      parameter (lvec = 600)
c
      implicit double precision (a-h,o-z)
c
      include 'common.h'
c
      dimension xl(8), vl(8), tl(4) 
c
      common /temp1 / ax(lvec)   , uaven(lvec),
     &                vaven(lvec), rhoan(lvec),
     &                rmfd(lvec) , rhoaven(lvec)
c
c..... array size check
c
      if (nelz+1 .gt. lvec)  then
        write(iout,1001)
        return
      end if
c
c..... zero arrays
c
      do 50 n=1, nelz+1
        ax(n) = zero
        rhoaven(n) = zero
        uaven(n) = zero
        vaven(n) = zero
50    continue
c
      if(abs(x(2, ien(1,1))-x(2, ien(2,1))).lt.eps7)then
c
c..... assuming elements increment across first
c
c..... localize coordinates, temperature and  velocity
c
      do 200 ix = 1, nelx 
c
c..... get the ix element for the iz row now
c
        do 100 iz = 1, nelz
          iel = ix + nelx*(iz - 1) 
          xl(1) = x(1,ien(iel,1) ) 
          xl(2) = x(1,ien(iel,2) ) 
c
          tl(1) = dens( ien(iel,1) ) 
          tl(2) = dens( ien(iel,2) ) 
          vl(1) = v(1, ien(iel,1) ) 
          vl(2) = v(2, ien(iel,1) ) 
          vl(3) = v(1, ien(iel,2) ) 
          vl(4) = v(2, ien(iel,2) ) 
          deltx = pt5 * ( xl(2) - xl(1) )
          ax(iz)    = ax(iz) + deltx * two
          rhoaven(iz) = rhoaven(iz) + (tl(1) + tl(2))*deltx
          uaven(iz) = uaven(iz)+(abs(vl(1))+abs(vl(3)))*deltx
          vaven(iz) = vaven(iz)+(abs(vl(2))+abs(vl(4)))*deltx
100     continue 
c
c Put the following inside the 200 loop and moved vl(nelz,5-8) and 
c tl(nelz,3-4) from 100 loop to save computation 6/1/94 jjita
c
        ax(nelz+1)    = ax(nelz) 
        tl(3) = dens(  ien(iel,3) ) 
        tl(4) = dens(  ien(iel,4) ) 
        rhoaven(nelz+1) = rhoaven(nelz+1) + 
     &                  (tl(3) + tl(4))*deltx
        vl(5) = v(1, ien(iel,3) ) 
        vl(6) = v(2, ien(iel,3) ) 
        vl(7) = v(1, ien(iel,4) ) 
        vl(8) = v(2, ien(iel,4) ) 
        uaven(nelz+1) = uaven(nelz+1) + 
     &                  (abs(vl(5)) + abs(vl(7)))*deltx
        vaven(nelz+1) = vaven(nelz+1) + 
     &                  (abs(vl(6)) + abs(vl(8)))*deltx
200   continue 
      else
c
c..... assuming elements increment vertical fastest
c
c..... localize coordinates, temperature and  velocity
c
      do 400 ix = 1, nelx 
c
c..... get the ix element for the iz row now
c
        do 300 iz = 1, nelz
          iel = iz + nelz*(ix - 1) 
          xl(1) = x(1,ien(iel,1) ) 
          xl(2) = x(1,ien(iel,2) ) 
c
          tl(1) = dens(  ien(iel,1) ) 
          tl(2) = dens(  ien(iel,2) ) 
          vl(1) = v(1, ien(iel,1) ) 
          vl(2) = v(2, ien(iel,1) ) 
          vl(3) = v(1, ien(iel,2) ) 
          vl(4) = v(2, ien(iel,2) ) 
          deltx = pt5 * ( xl(2) - xl(1) )
          ax(iz)    = ax(iz) + deltx * two
          rhoaven(iz) = rhoaven(iz) + (tl(1) + tl(2))*deltx
          uaven(iz) = uaven(iz)+(abs(vl(1))+abs(vl(3)))*deltx
          vaven(iz) = vaven(iz)+(abs(vl(2))+abs(vl(4)))*deltx
300     continue 
c
c Put the following inside the 400 loop and moved vl(nelz,5-8) and
c tl(nelz,3-4) from 300 loop to save computation 6/1/94 jjita
c
        ax(nelz+1)    = ax(nelz)
        tl(3) = dens(  ien(iel,3) ) 
        tl(4) = dens(  ien(iel,4) ) 
        rhoaven(nelz+1) = rhoaven(nelz+1) + 
     &                  (tl(3) + tl(4))*deltx
        vl(5) = v(1, ien(iel,3) ) 
        vl(6) = v(2, ien(iel,3) ) 
        vl(7) = v(1, ien(iel,4) ) 
        vl(8) = v(2, ien(iel,4) ) 
        uaven(nelz+1) = uaven(nelz+1) + 
     &                  (abs(vl(5)) + abs(vl(7)))*deltx
        vaven(nelz+1) = vaven(nelz+1) + 
     &                  (abs(vl(6)) + abs(vl(8)))*deltx
400     continue 
      end if
c
c..... divide mean quantities by the row length
c
        do 500 iz = 1, nelz 
         rhoaven(iz)   =  rhoaven(iz)/ax(iz)
         uaven(iz)   =  uaven(iz)/ax(iz)
         vaven(iz)   =  vaven(iz)/ax(iz)
500     continue
        rhoaven(nelz+1)   =  rhoaven(nelz+1)/ax(nelz+1)
        uaven(nelz+1)     =  uaven(nelz+1)/ax(nelz+1)
        vaven(nelz+1)     =  vaven(nelz+1)/ax(nelz+1)
c
c...calculate radial mass flux diagnostic
c
      write(icomp,*) "radial mass flux", lstep, time
      scale = 0.0e0
      xscale = 0.0e0
      do 900 iz = 1, nelz
        iel = iz
        xscale = xscale + shl(1,5)*rhoaven(iz)*vaven(iz)*det(iel,1)
     &              + shl(4,5)*rhoaven(iz+1)*vaven(iz+1)*det(iel,4)
        scale = scale + shl(1,5)*det(iel,1) + shl(4,5)*det(iel,4)
900   continue 
c
c Inverted scale to save division in this loop on 6/1/94 jjita
c
      scale = scale / xscale
      do 950 iz = 1, nelz
      if(abs(x(2,ien(1,1))-x(2,ien(2,1))).lt.eps7) then
          iel = 1 + (iz-1)*nelx 
      else
          iel = iz
      endif
          xl(2) = x(2,ien(iel,1) ) 
          rmfd(iz) = rhoaven(iz)*vaven(iz) * scale
          write(icomp,1000) xl(2),rmfd(iz)
950   continue
      write(icomp,1000) one,rmfd(1)
c
c..... return
c
 1000 format(f6.4,2x,f6.4)
1001  format( ' WARNING: Dimension problem in subroutine mean - ',//,
     &        '  -->  no mean calculations performed' )
      return
      end


