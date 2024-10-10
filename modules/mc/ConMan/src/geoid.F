      subroutine geoid (x, t, stress)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                  c
c   subroutine to calculate GEOID from buoyancy and stress field   c
c   from ConMan in the 2-D Cartesian case. GEOID uses the          c
c   runfile to get the geom file (for coordinates) the             c
c   temperature & velocity file for the buoyancies and the         c
c   stress file to get the topography at the surface and bottom.   c
c   Written by Scott King Feb 9, 1990 at Caltech.                  c
c                                                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   
c   parameter for array dimensions, greater than the number of nodes
c 
c
      implicit double precision (a-h,o-z)
      PARAMETER (NWAVE=2000)
c
#include<common.h>
#include<materials.h>
c
      double precision x(2,*), t(*),  stress(6,*)
      real    z1, z2, z3, aratio, pi
      real    tlay(NWAVE), tintz(NWAVE), wave(NWAVE,NWAVE)
      real    geoid1(NWAVE),rtzzt(NWAVE),rtzzb(NWAVE), geoid2(NWAVE)
c
c constants for dimensionalization (geoid and topography in meters)
c
      pi = 3.14159272e0
      aratio = x(1,numnp)
c
c  The values used for the Blankenback et al. Benchmark 1a, 1b, 1c
c  remove the comments from the following eight lines of code and
c  comment out the lines below for the Blankenbach benchmarks
c
c     d = 1.0e6
c     g = 1.0e1
c     rho = 4.0e3
c     evisc = 2.5e19 * rho * 1.0e4 / ra(1)
c     dif = 1.0e-6
c     bigG = 6.673e-11
c     alfa = 2.5e-5
c     deltaT = 1.0e3
c
c  The values below are more representative of whole mantle convection
c  values.  Notice that this combination of values assumes a specific
c  values of the Rayleigh number.   If the input value of the Rayleigh
c  number differs, then these values will not give meaningful geoid
c  and topography values.
c
c  These lines should be commented out for the Blankenbach problems.
c
      d = 3.000e6
      g = 1.0e1
      rho = 4.0e3
      evisc = 4.32e21
      dif = 1.0e-6
      bigG = 6.673e-11
      alfa = 2.0e-5
      deltaT = 2.000e3
c
c  Below here nothing needs to be changed.
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
      if (nwrap .gt. 0) npoints = nelx 
      if (nwrap .gt. 0) ntrans = nelx 
      if (nwrap .eq. 0) then 
          do 101 i=1,npoints
           rtzzt(i) =  -stress(5,(nelz+1)*i)
           rtzzb(i) =  -stress(5,1+(nelz+1)*(i-1))
           rtzzt(2*npoints-i) =  rtzzt(i)
           rtzzb(2*npoints-i) =  rtzzb(i)
101       continue
          ntrans = 2*npoints - 2
      endif
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
          do 600 iz = 1, nelz+1
            do 500 ix = 1, npoints
              tlay(ix) = t(iz + (nelz+1) * (ix-1))
500         continue
            if (nwrap .eq. 0) then 
              do 501 i=1,npoints
                tlay(2*npoints-i) = tlay(i)
501           continue
            endif
            call fftl(tlay, ntrans, 1, ierror)
c
c   move twave
c
            do 550 ik=1,ntrans
              wave(iz,ik) = tlay(ik)
550         continue
600       continue
c                         -kz
c   integrate with depth e 
c   this assumes a z fastest node ordering!
c
          numwav = ntrans/2 + 1
          do 900 ik = 1, numwav
c
c integral over shape by Simpson's Rule
c
c note that "jk" is the array index while "ik" is the wavenumber
c (there are two coeficients for each wavenumber
c
            jk = 2*ik
            tintz(jk-1) = 0.0e0
            tintz(jk)   = 0.0e0
c
c use a 3 point 3rd order rule for the integral
c
c          do 800 iz = 2, nelz, 2
            do 800 iz = 2, nelz+1
              z3 = x(2,iz+1)
              z2 = x(2,iz)
              z1 = x(2,iz-1)
c             tintz(jk-1) = tintz(jk-1) + 
c    &   ( exp(-1.0*pi*float(ik-1)*(1.0-z3)/aratio)*wave(iz+1,jk-1)+
c    & 4.0*exp(-1.0*pi*float(ik-1)*(1.0-z2)/aratio)*wave(iz,jk-1)+
c    &     exp(-1.0*pi*float(ik-1)*(1.0-z1)/aratio)*wave(iz-1,jk-1))*
c    &      (z2-z1)/3.0
c             tintz(jk) = tintz(jk) + 
c    &   ( exp(-1.0*pi*float(ik-1)*(1.0-z3)/aratio)*wave(iz+1,jk)+
c    & 4.0*exp(-1.0*pi*float(ik-1)*(1.0-z2)/aratio)*wave(iz,jk)+
c    &     exp(-1.0*pi*float(ik-1)*(1.0-z1)/aratio)*wave(iz-1,jk))*
c    &      (z2-z1)/3.0
c two point rule
              tintz(jk-1) = tintz(jk-1) +
     &  ( exp(-1.0*pi*float(ik-1)*(one-z2)/aratio)*wave(iz,jk-1)+
     &    exp(-1.0*pi*float(ik-1)*(one-z1)/aratio)*wave(iz-1,jk-1))*
     &  (z2-z1)/2.0
              tintz(jk) = tintz(jk) +
     &  ( exp(-1.0*pi*float(ik-1)*(one-z2)/aratio)*wave(iz,jk)+
     &    exp(-1.0*pi*float(ik-1)*(one-z1)/aratio)*wave(iz-1,jk))*
     &  (z2-z1)/2.0
800         continue
c
900       continue
c
c zero out the k=0 harminics to remove the mean
c
        rtzzt(1)  = 0.0
        rtzzt(2)  = 0.0
        rtzzb(1)  = 0.0
        rtzzb(2)  = 0.0
        tintz(1)  = 0.0
        tintz(2)  = 0.0
        geoid1(1) = 0.0
        geoid1(2) = 0.0
c
c   spit out the result
c
        do 1100 ik = 2, numwav
        jk = 2*ik
        geoid1(jk-1) = (1.0/float(ik-1))*(rtzzt(jk-1)*evisc*dif/(d*d*g) 
     &  -rtzzb(jk-1)*exp(-1.0*pi*float(ik-1)/aratio)*evisc*dif/(d*d*g)
     &  -tintz(jk-1)*rho*alfa*deltaT*d )
        geoid1(jk)   = (1.0/float(ik-1))*(rtzzt(jk)*evisc*dif/(d*d*g) 
     &  -rtzzb(jk)*exp(-1.0*pi*float(ik-1)/aratio)*evisc*dif/(d*d*g)
     &  -tintz(jk)*rho*alfa*deltaT*d )
        geoid2(jk-1)=-(1.0/float(ik-1))*
     &  tintz(jk-1)*rho*alfa*deltaT*d
        geoid2(jk)  =-(1.0/float(ik-1))*
     &  tintz(jk)*rho*afha*deltaT*d
1100    continue
c
1006  format(1x,i5,1x,2(1pe9.3,1x),5(1pe12.5,1x))
1007  format(1x,i5,1x,5(1pe12.5,1x))

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Transforms Wave Number Data Back to real world for use with plot
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call fftl(rtzzt,  ntrans, 2, ierror)
      call fftl(rtzzb,  ntrans, 2, ierror)
      call fftl(geoid1, ntrans, 2, ierror)
      call fftl(geoid2, ntrans, 2, ierror)
c
c warning if npoints is odd, the last point is garbage
c
      pk = 2 * bigG*d/g
      do 200 i=1,npoints
        xm = aratio*float(i-1)/(npoints-1)
        write(igeoid,1001) xm, rtzzt(i)*evisc*dif / (rho*d*d*g) ,
     &  pk*geoid2(i), pk * geoid1(i) 
200   continue
c     call myflush(igeoid)
c
1001  format(4e15.7)
      return
      end
