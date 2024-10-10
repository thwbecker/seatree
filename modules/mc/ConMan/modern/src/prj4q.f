      subroutine prj4q  (xl,     yl,     gxl,    gyl,  soll,   npgrd,  
     &                   psoll,  lpsoll, iel)
c
c-----------------------------------------------------------------------
c
c  This routine compute the shapefunctions for a 4-node bilinear
c quadrilateral at a set of given global points.  Then projects the 
c solution of the element into those points.
c
c input:
c  xl:     X value of the nodes
c  yl:     Y value of the nodes
c  gxl:    X where the shape-function to be computed
c  gyl:    Y where the shape-function to be computed
c  soll:   solution at the nodes
c  npgrd:  number of points to be projected to
c  ncol:   number of data columns
c  iel:    element number (for error)
c
c output:
c  psoll: solution at the projected points
c  lpsoll: flag indicating the point in element
c
c
c algorithm:
cconstruct the system:
c  x = ax + bx xi + cx eta + dx xi eta
c  y = ay + by xi + cy eta + dy xi eta
c where ax, bx, cx, dx, ay, by, cy and dy are the element
c geometrical constants and xi and eta are the local coordinates.  
c Solve this system for xi and eta for a given x and y, 
c from which compute the shape-functions.
c
c Note: a point is considered inside a bilinear element 
c  f xi=[-tol,tol] and eta=[-tol,tol], where tol=1+eps.
c
c Number of operations:
c  28  flop per element
c  50  flop per projection point
c	 7  flop per projection point per data column
c
c Farzin Shakib, Spring 1989.
c-----------------------------------------------------------------------
c
        implicit double precision (a-h,o-z)
        parameter (MAXPNT = 8000000)
        include 'common.h'
c
        dimension    xl(*),    yl(*),    gxl(*),   gyl(*),   
     &               soll(*),     psoll(*),    lpsoll(*)
c
        dimension    xi(MAXPNT),      eta(MAXPNT),     a11(MAXPNT),
     &               a12(MAXPNT),     a21(MAXPNT),     a22(MAXPNT),
     &               dinv(MAXPNT),    xip(MAXPNT),     xim(MAXPNT),
     &               etap(MAXPNT),    etam(MAXPNT),    shp(4,MAXPNT)
c
c.... set the tolerance
c
      tol = one + 1.0d-3
c
c.... construct the system
c
      ax = (xl(1) + xl(2)) + (xl(3) + xl(4))
      bx = (xl(2) + xl(3)) - (xl(1) + xl(4))
      cx = (xl(3) + xl(4)) - (xl(1) + xl(2))
      dx = (xl(1) + xl(3)) - (xl(2) + xl(4))
c
      ay = (yl(1) + yl(2)) + (yl(3) + yl(4))
      by = (yl(2) + yl(3)) - (yl(1) + yl(4))
      cy = (yl(3) + yl(4)) - (yl(1) + yl(2))
      dy = (yl(1) + yl(3)) - (yl(2) + yl(4))
c
c.... eliminate explicit dependence of system to ax and ay
c
      do k = 1, npgrd
        gxl(k) = four * gxl(k) - ax
        gyl(k) = four * gyl(k) - ay
      enddo
c
c.... get a linear solution (i.e., ignore the xi eta term)
c
      detinv = bx * cy - cx * by
      if (detinv .lt. eps7) goto 999
      detinv = one / detinv
c
      do k = 1, npgrd
        xi(k)  = detinv * (cy * gxl(k) - cx * gyl(k))
        eta(k) = detinv * (bx * gyl(k) - by * gxl(k))
      enddo
c
c.... do one Newton correction
c
      itmp = 0
      do 300 k = 1, npgrd
        gxl(k)  = gxl(k) + dx * (xi(k) * eta(k))
        gyl(k)  = gyl(k) + dy * (xi(k) * eta(k))
c
        a11(k)  = bx + dx * eta(k)
        a12(k)  = cx + dx * xi(k)
        a21(k)  = by + dy * eta(k)
        a22(k)  = cy + dy * xi(k)
c
        dinv(k) = (a11(k) * a22(k)) - (a12(k) * a21(k))
        if (dinv(k) .lt. eps7) itmp = 1
 300  continue
c
       if (itmp .eq. 1) goto 999
c
       do 400 k = 1, npgrd
         dinv(k) = one / dinv(k)
c
         xi(k)   = dinv(k) * (a22(k) * gxl(k) - a12(k) * gyl(k))
         eta(k)  = dinv(k) * (a11(k) * gyl(k) - a21(k) * gxl(k))
400    continue
c
c.... compute the shape functions
c
      do 500 k = 1, npgrd
        xip(k)  = pt5 * (one + xi(k))
        xim(k)  = pt5 * (one - xi(k))
        etap(k) = pt5 * (one + eta(k))
        etam(k) = pt5 * (one - eta(k))
c
        shp(1,k) = xim(k) * etam(k)
        shp(2,k) = xip(k) * etam(k)
        shp(3,k) = xip(k) * etap(k)
        shp(4,k) = xim(k) * etap(k)
500   continue
c
c.... project the solution
c
      do 600 k = 1, npgrd
        psoll(k) = shp(1,k) * soll(1) + shp(2,k) * soll(2) +
     &    shp(3,k) * soll(3) + shp(4,k) * soll(4)
600   continue
c
c.... set the projection flag (if point is inside the element)
c
      do 800 k = 1, npgrd
        lpsoll(k) = 0
        if (dabs(xi(k)).le.tol .and. dabs(eta(k)).le.tol) lpsoll(k) = 1
800   continue
c
c.... return
c
      return
c
c.... error handling
c
 999  write(6,*) "prj4q: zero determanent of Jacobian; element:", iel
c
c.... end
c
      end
