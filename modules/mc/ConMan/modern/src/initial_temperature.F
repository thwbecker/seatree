      subroutine  intemp( x ,  t , tdot)
c
      use bc
      implicit double precision (a-h,o-z)
c
      include 'common.h'
c
      dimension   x(nsd,*) , t(*) , tdot(*)

      pix = 3.14159265358979323D+0/xsize
      piz = 3.14159265358979323D+0/zsize

c.... coefficient for error function
c.... smaller zcoeff means a thick initial stagnat lid
      zcoeff = 5.00000000000000D+0

      if ( inrstr .eq. 0 ) then
c
c.... use perturbed conductive start
c

        do 100 i = 1,numnp
          tdot(i) = zero
          k = idt(i)
100     continue

c.... maximum of the bottom temperature (when z equals 0)
        
        bmax = one - T0 / diff_T * (exp(Di) - one)
        do 200 i = 1,numnp
          k = idt(i)
          if(k .ne. 0) then
c
c....[1]error function temperature profile (top 0, bottom 1)
c          zvalue = (one - x(2,i)) * zcoeff
c          call efunc(zvalue,tvalue)

cc....  combination of sin and cos functions
c          t(i) = tvalue * bmax
c     &           + pert * dcos(x(1,i)*pix) * dsin(x(2,i)*piz)        
c..............................................................

c....[2]linear temperature profile (top 0, bottom 1)

cc....  combination of sin and cos functions
c           t(i) = (zsize-x(2,i)) / zsize * bmax
            t(i) = 0.60d0 * bmax
     &           + pert * dcos(x(1,i)*pix) * dsin(x(2,i)*piz)        
c       tzmin = 0.0
c       tzbl = 0.05
c       tPI=4.0*ATAN(1.0)
c       tT0 = 0.6
c       tamplitude = 0.05
c         tz = (1.0-tt0)*DEXP(-(1.0-(x(2,i)-tzmin))/tzbl)
c    &             -tt0*(DEXP(-(x(2,i)-tzmin)/tzbl)-1.0)
c         tpert= tamplitude*DCOS(2.0*tPI*x(1,i))*DSIN(2.0*tPI*x(2,i))
c         t(i) = MAX(MIN(1.0,tz+tpert),0.0)
c
c....[3]adiabatic temperature profile
c            t(i) = T0 / diff_T * (exp(Di*(one-x(2,i))) - one)
c     &           + pert * dcos(x(1,i)*pix) * dsin(x(2,i)*piz)

          else
            t(i) = tbc(i)
          endif
        
200     continue
       
       time0=zero
       istep0=0
      else
c
c.... use restart file
c
        read(irsin,1000) istep0, time0
        read(irsin,*)
          do 301 n = 1, numnp
            read(irsin,1001) k,t(n),tdot(n)  
301       continue
      end if
1000  format(7x,i7,10x,f10.6)
1001  format(1x,i7,4(1pe15.8,1x))
c
c.... return
c
      return    
      end
