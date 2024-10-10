      subroutine geoid (x, t, tzzt, tzzb)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                  c
c   subroutine to calculate GEOID from buoyancy and stress field   c
c   from ConMan in the 2-D Cartesian case. GEOID uses the          c
c   runfile to get the geom file (for coordinates) the             c
c   temperature & velocity file for the buoyancies and the         c
c   stress file to get the topography at the surface and bottom.   c
c   Written by Scott King Feb 9, 1990 at Caltech.                  c
c                                                                  c
c   density variation added and bugs are fixed.                    c
c   Revised by Changyeol Lee May, 2009 at Virginia Tech.           c
c                                                                  c
c   non-dimensionalized input values                               c
c   x(2, numnp)    : x and z coordinates                           c
c   t(numnp)       : temperature                                   c
c   tzzt(nelx + 1) : Tau zz of dynamic stress at top of the box    c
c   tzzb(nelx + 1) : Tau zz of dynamic stress at bottom of the box c
c   constant for dimensionalization (geoid and topography in meters)
c   temperature & pressure dependent viscosity is already 
c   considered in the subroutine of stress.F as a evisc(intp)
c   d      : depth of the box (m)
c   g      : gravity acceleration (m/s2)
c   rho    : reference density (Kg/m3) (multiplied by density)
c   visc0  : reference viscosity (Pa.S)
c   diff   : thermal diffusivity (m2/s)
c   bigG   : universial constant of gravity (N.m2/kg2)
c   expans : thermal expansivity (/K)
c   diff_T : temperature difference between the top and the bottom (K)
c   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      use materials, only : tref, dens, ra
      implicit double precision (a-h,o-z)
c   
c   parameter for array dimensions, greater than the number of nodes
      parameter (NWAVE = 10000)
c 
      include 'common.h'
c
      real z1, z2, z3, aratio, pi
      double precision x(2,*), t(*), tzzt(*), tzzb(*)
      real tlay(NWAVE), tintz(NWAVE), wave(NWAVE,NWAVE)
      real geoid1(NWAVE), rtzzt(NWAVE), rtzzb(NWAVE), geoid2(NWAVE),
     &     geoid3(NWAVE), geoid4(NWAVE) 
c
c
c constants for dimensionalization (geoid and topography in meters)
c  the values used for the Blankenback et al. Benchmark 1a
c   (Ra=1.0e4, constant viscosity)
c 
c......... set for Suite2 for TALA benchmarking 
c      d      = 3.00e+6                  ! 3000 km
c      g      = 1.00e+1                  
c      rho    = 4.00e+3
c      visc0   = 1.00e+23
c      diff   = 1.00e-6
c      expans = 3.0e-5
c      expans = 6.0e-5    ! case2a
c      expans = 9.0e-5    ! case2b
c      expans = 1.2e-4    ! case2c
c.....................................    

c......... set for Blankenbach et al.
      d      = 1.00e+6          
      g      = 1.00e+1                  
      visc0  = 1.00e+27/ra(1)          
      diff   = 1.00e-6
      expans = 2.5e-5
      rho    = 4.00e+3
      rhotop = rho    ! density contrast at top
      rhobot = rho    ! density contrast at bottom

c.....................................    

      bigG   = 6.673e-11
      pi = 4.0*atan(1.0)
      aratio = 2.0e0 * x(1,numnp)        ! wavelength

      if (nwrap .gt. 0) then
       aratio = x(1,numnp)               ! wavelength for wrapped BC
      endif 
c
c.... initialization
c
      do i = 1, NWAVE 
        tlay(i) = 0.0e0  ; tintz(i) = 0.0e0 ; geoid1(i) = 0.0e0
        rtzzt(i) = 0.0e0 ; rtzzb(i) = 0.0e0 ; geoid2(i) = 0.0e0
        geoid3(i) = 0.0e0 ; geoid4(i) = 0.0e0
        do j = 1, NWAVE
          wave(i,j) = 0.0e0
        enddo
      enddo
c
c if the problem uses reflecting boundary conditions, this is correct
c  we also need to reflect the solution about the right boundary
c  to get the correct symmetry.
c
      npoints = nelx + 1
c
c otherwise discard the last point - it is really the first (not tested)
c   nb: if npoints is odd, the tranform will ignore the last point
c 
      ntrans = 2*npoints - 2

      if (nwrap .gt. 0) npoints = nelx
      if (nwrap .gt. 0) ntrans = nelx 
      do i=1,npoints
         rtzzt(i) = tzzt(i)
         rtzzb(i) = tzzb(i)
         if (nwrap .eq. 0) then 
           rtzzt(2*npoints-i) = tzzt(i)
           rtzzb(2*npoints-i) = tzzb(i)
         endif
      enddo
c
c   change sign of rtzzt and tzzt
c   at the top : compression (stress <0) -> positive topography : change sign
c   at the bottom : compression (stress <0) -> negative topography : leave sign
      do i = 1, (2*npoints - 1)       
         rtzzt(i) = - rtzzt(i)
         rtzzb(i) = rtzzb(i)
      enddo  
      do i = 1, npoints
         tzzt(i)  = - tzzt(i)
         tzzb(i)  =  tzzb(i)
      enddo
c
c   calculate fourier decomposition of the top and bottom
c   nb: if npoints is odd, the tranform will ignore the last point
c
      call fftl(rtzzt, ntrans, 1, ierror)
c
c   bottom: (assumes iz is the fastest varying index)
c
      call fftl(rtzzb, ntrans, 1, ierror)
c
c   calculate fourier decomposition at each layer
c
      do iz = 1, nelz+1
        do ix = 1, npoints
          node  = iz + (nelz+1)*(ix-1)
          tlay(ix) = t(node) * dens(node) + tref(node) 
        enddo
        if (nwrap .eq. 0) then 
          do i=1,npoints
            tlay(2*npoints-i) = tlay(i)
          enddo
        endif
        call fftl(tlay, ntrans, 1, ierror)
c
c   move twave
c
        do ik=1,ntrans
           wave(iz,ik) = tlay(ik)
        enddo 
      enddo
c                         -kz
c   integrate with depth e 
c   this assumes a z fastest node ordering!
c   integral over shape by Simpson's Rule
c
      numwav = ntrans/2 + 1
      do ik = 1, numwav
c note that "jk" is the array index while "ik" is the wavenumber
c (there are two coeficients for each wavenumber
        jk = 2*ik
        tintz(jk-1) = zero
        tintz(jk)   = zero
        do iz = 2, nelz+1
          z3 = x(2,iz+1)
          z2 = x(2,iz)
          z1 = x(2,iz-1)
c two point rule
          tintz(jk-1) = tintz(jk-1) +
     &                ( exp(-2.0e+0*pi*(ik-1)*(one-z2)/aratio)
     &                 *wave(iz,jk-1)+
     &                  exp(-2.0e+0*pi*(ik-1)*(one-z1)/aratio)
     &                 *wave(iz-1,jk-1))
     &                 *(z2-z1)/2.0e+0
          tintz(jk)   = tintz(jk) +
     &                ( exp(-2.0e+0*pi*(ik-1)*(one-z2)/aratio)
     &                 *wave(iz,jk)+
     &                  exp(-2.0e+0*pi*(ik-1)*(one-z1)/aratio)
     &                 *wave(iz-1,jk))
     &                 *(z2-z1)/2.0e+0
        enddo
      enddo
c
c zero out the k=0 harmonics to remove the mean
c
      rtzzt(1)  = zero
      rtzzt(2)  = zero
      rtzzb(1)  = zero
      rtzzb(2)  = zero
      tintz(1)  = zero
      tintz(2)  = zero
      geoid1(1) = zero
      geoid1(2) = zero
      geoid2(1) = zero
      geoid2(2) = zero
      geoid3(1) = zero
      geoid3(2) = zero
      geoid4(1) = zero
      geoid4(2) = zero
c
c   spit out the result
c
      do ik = 2, numwav
        jk = 2*ik
c    visc0*diff/(d*d*g*rhotop) : height of dynamic topography, top
c    visc0*diff/(d*d*g*rhobot) : height of dynamic topography, bottom
c... whole contribution
        geoid1(jk-1) = (1.0e0/(ik-1))* 
     &  ( rtzzt(jk-1)*visc0*diff/(d*d*g*rhotop)*rhotop ! contribution of excess mass above the top
     &   +rtzzb(jk-1)*exp(-2.0e0*pi*(ik-1)/aratio)*    ! contribution of excess mass under the bottom
     &   visc0*diff/(d*d*g*rhobot)*rhobot
     &   -tintz(jk-1)*rho*expans*diff_T*d)             ! contribution of deficit mass due to thermal expansion
        geoid1(jk)   = (1.0e0/(ik-1))*
     &  ( rtzzt(jk)*visc0*diff/(d*d*g*rhotop)*rhotop
     &   +rtzzb(jk)*exp(-2.0e0*pi*(ik-1)/aratio)*
     &   visc0*diff/(d*d*g*rhobot)*rhobot
     &   -tintz(jk)*rho*expans*diff_T*d) 
c... temperature contribution
        geoid2(jk-1)=-(1.0e0/(ik-1))*
     &  tintz(jk-1)*rho*expans*diff_T*d 
        geoid2(jk)  =-(1.0e0/(ik-1))*
     &  tintz(jk)*rho*expans*diff_T*d 
c... bottom contribution
        geoid3(jk-1)=(1.0e0/(ik-1))*rtzzb(jk-1)
     &  *exp(-2.0e0*pi*(ik-1)/aratio)*visc0*diff/(d*d*g*rhobot)*rhobot
        geoid3(jk)  =(1.0e0/(ik-1))*rtzzb(jk)
     &  *exp(-2.0e0*pi*(ik-1)/aratio)*visc0*diff/(d*d*g*rhobot)*rhobot
c... top contribution
        geoid4(jk-1)=(1.0e0/(ik-1))*rtzzt(jk-1)
     &              *visc0*diff/(d*d*g*rhotop)*rhotop
        geoid4(jk)  =(1.0e0/(ik-1))*rtzzt(jk)
     &              *visc0*diff/(d*d*g*rhotop)*rhotop

      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Transforms Wave Number Data Back to real world for use with plot
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call fftl(rtzzt,  ntrans, 2, ierror)
      call fftl(rtzzb,  ntrans, 2, ierror)
      call fftl(geoid1, ntrans, 2, ierror)
      call fftl(geoid2, ntrans, 2, ierror)
      call fftl(geoid3, ntrans, 2, ierror)
      call fftl(geoid4, ntrans, 2, ierror)
c
c warning if npoints is odd, the last point is garbage
c
      write(igeoid,1002)
       
      pk = 2.0d0 * pi * bigG / (2.0d0 * pi / aratio * g) * d

      xlength = x(1,numnp)
      if (nwrap .gt. 0) then
      xlength = x(1,numnp) - x(1,numnp)/nelx  ! wrapped BC
      endif

      do 200 i=1,npoints
        xm = xlength*float(i-1) / (npoints-1)
        write(igeoid,1001) xm, time,
     &  rtzzt(i)*visc0*diff / (rhotop*d*d*g) ,
     &  rtzzb(i)*visc0*diff / (rhobot*d*d*g) , 
     &  pk * geoid1(i) 
200   continue
1001  format(5e15.7)
1002  format(' x coordinate        time     top topo       bottom topo      
     &      geoid  ')
      return
      end
