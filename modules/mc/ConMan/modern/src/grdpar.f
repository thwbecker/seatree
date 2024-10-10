      subroutine grdpar (x, gpars,  ngpars)
c
c-----------------------------------------------------------------------
c
c  This routine reads the grid parameters.
c
c
c input:
c  x            : unstructured mesh coordinates
c
c output:
c  gpars      : structured mesh parameters (xmin, xmax, xinc,
c                                    ymin, ymax, yinc)
c  ngpars      : structured mesh resolution parameters 
c
c Farzin Shakib, Spring 1989.
c-----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
c     maximum pixel for x and z direction is 301 
      parameter (ipsolx = 301)
      include 'common.h'
c  
      dimension x(2,*), gpars(3,2), ngpars(2)
c
c.... compute the maximum and minimum mesh values
c
      xmin = x(1,1)
      xmax = x(1,numnp)
      zmin = x(2,1)
      zmax = x(2,numnp)
c
c.... get grid resolution
c
      ngpnt = ipsolx - 1
      fact  = dmin1(one, (xmax-xmin)/(zmax-zmin))
      ngpnt = int(fact * dble(ngpnt))
      ngpars(1)  = ngpnt + 1
      hxz        =     (xmax - xmin) / dble(ngpnt)
      itmp       = int((zmax - zmin) / hxz + 0.5)
      ngpars(2)  = itmp
      gpars(1,1) = xmin
      gpars(2,1) = xmax
      gpars(1,2) = zmin
      gpars(2,2) = zmax
      gpars(3,1) = hxz
      gpars(3,2) = hxz
      gpars(2,2) = gpars(1,2) + float(itmp) * hxz
c
c.... return
c
      return
      end
