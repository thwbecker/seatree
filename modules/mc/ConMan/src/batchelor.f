      subroutine batchelor(xorig, yorig, u, v)
c
c  A program to implement the Batchelor cornerflow solution
c  See Batchelor (1967) pp. 224-225
c
c  Written: S. D. King November 1, 2002
c
c  constants
c
      implicit double precision (a-h,o-z)
      common /io    / iin,    igeom, iout , itsout , itout , imout,
     &                irsin , irsout

      two = 2.0e0
      four = 4.0e0
      pi = four*atan(1.0)
c
c velocity of the downgoing slab
c
c     Uplate  = 1.538330
      Uplate  = 1.54
      Uplate  = 1.53
      visc = 1.0e3
c
c the general expression for the streamfunction is given by:
c    psi = Ax + By + (Cx + Dy)*arctan(y/x)
c
c and 
c    u   = - d psi/dy = -B -x*(Cx + Dy)/(x^2+y^2) - D arctan(y/x)
c    
c    v   =   d psi/dx = A - y*(Cx + Dy)/(x^2+y^2) + C arctan(y/x)
c
c for the arc corner with a slab at theta = pi/4 where the downgoing 
c  slab velocity is U*sqrt(two)/two, we have
c
      A = 0.0e0
      B =  pi*Uplate/(two-(pi**2)/four)
      C = -pi*Uplate/(two-(pi**2)/four)
      D = -Uplate*two*(two-pi/two)/(two-(pi**2)/four)
c
c   works for my grids
      x = 1.0e3*(610.0e0-xorig)
      y = 1.0e3*(550.0e0-yorig) 
c   works for changyeol's grids
c     x = 1.0e3*(xorig-50.0e0)
c     y = 1.0e3*(550.0e0-yorig) 
      press = -2*visc*(C*x + D*y)/((x**2)+(y**2))
      psi   = A*x + B*y + (C*x + D*y)*atan2(y,x)
      u     = -B - x*(C*x + D*y)/((x**2)+(y**2)) 
     &           - D*atan2(y,x)
      v     =  A - y*(C*x + D*y)/((x**2)+(y**2)) 
     &           + C*atan2(y,x)
c     if (u .le. Uplate) u = Uplate
c     if (v .le. Uplate) v = Uplate
c needed for sdk grids
      u = -u
      v = -v
      if (y .le. 0) then
        u = 0.0d0
        v = 0.0d0
      endif
c     write(iout,*) xorig,yorig,x,y,u,v
c
c end
c
      return
      end
