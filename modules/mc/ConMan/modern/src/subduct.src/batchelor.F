c
c     box height (600)
c     plate thickness (50)
c     plate velocity (2.1637467D+0)
c      
      subroutine batchelor(xorig, yorig, u, v, box_height,
     &     plate_thickness, Uplate)
c
c  A program to implement the Batchelor cornerflow solution
c  See Batchelor (1967) pp. 224-225
c
c  Written: S. D. King November 1, 2002
c
c  constants
c
      implicit double precision (a-h,o-z)
c
      double precision box_height, plate_thickness, Uplate, poff
      
      two = 2.0D+0
      four = 4.0D+0
      pi = four*atan(1.0D+0)

C     offset, typically 600-50 = 550
      poff = box_height - plate_thickness;
c
c velocity of the downgoing slab
c Uplate is not x and y components
c      Uplate  = 2.1637467D+0 
      visc = 1.0D+3
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
      A = 0.0D+0
      B =  pi*Uplate/(two-(pi**2)/four)/sqrt(two)
      C = -pi*Uplate/(two-(pi**2)/four)/sqrt(two)
      D = -Uplate*sqrt(two)*(two-pi/two)/(two-(pi**2)/four)
c
c correct for undeformed grid
c     x = 1.0D+0*(xorig-50.000000D+0)
      x = 1.0D+0*(xorig - plate_thickness)
      
c correct for deformed grid
c     x = 1.0D+0*(610.0D+0-xorig)
      y = 1.0D+0*(poff - yorig)
      press = -2.0D+0*visc*(C*x + D*y)/((x**2)+(y**2))
      psi   = A*x + B*y + (C*x + D*y)*atan2(y,x)
      u     = -B - x*(C*x + D*y)/((x**2)+(y**2)) 
     &           - D*atan2(y,x)
      v     =  A - y*(C*x + D*y)/((x**2)+(y**2)) 
     &           + C*atan2(y,x)
c correct for undeformed grid
      u =  u
      v = -v
c correct for deformed grid
c     u = -u
c     v = -v
      if ((y .lt. 0).or.(x .lt. 0)) then
        u = 0.0d0
        v = 0.0d0
      endif
c
c don't set the top of the slab
c
      if (abs(x-y) .lt. 0.1) then
        u = 0.0d0
        v = 0.0d0
      endif
c     if ((abs(yorig-(xorig-60.0)) .lt. 1.0)) then
c       write(6,*) x,y,xorig,yorig,u,v
c     endif
c
c end
c
      return
      end
