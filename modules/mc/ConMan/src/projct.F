      subroutine projct (sol, psol ) 
c
c-----------------------------------------------------------------------
c
c  This routine projects the solution from an unstructured mesh to
c a structured grid.
c
c input:
c  sol		: solution on unstructured mesh
c  npsol	: array size of grid solution
c
c output:
c  psol		: projected solution
c  lpsol	: flag to indicate a projected solution
c
c Farzin Shakib, Spring 1989.
c-----------------------------------------------------------------------
c
       implicit double precision (a-h,o-z)
       parameter (MAXPNT = 111)
#include<common.h>
#include<fe.h>
#include<variables.h>
c
       dimension    sol(*),  psol(111,101),  lpsol(111,101)
c
       dimension    xl(4),  yl(4),  gxl(MAXPNT), gyl(MAXPNT),
     &            soll(4),  psoll(MAXPNT), lpsoll(MAXPNT)
c
       integer maxi, maxj
c
c.... setup some constants
c
      x0   = x(1,1)
      y0   = x(2,1)
c
      maxi = 111
      maxj = 101
c
c.... initialize the lpsol flag
c
      do i = 1, maxi
        do j = 1, maxj
          lpsol(i,j) = 0
        enddo
      enddo
c
c.... loop through the elements
c
      do 2000 iel = 1, numel
c
c.... compute the element parameters
c
        nenl = 4
c
        do n = 1, nenl
          node     = ien(iel,n)
          xl(n)    = x(1,node)
          yl(n)    = x(2,node)
          soll(n) = sol(node)
        enddo
c
c.... compute the rectangular block surrounding the element
c
        xmin = xl(1)
        xmax = xl(1)
        ymin = yl(1)
        ymax = yl(1)
c
        do n = 2, nenl
          xmin = min (xmin, xl(n))
          xmax = max (xmax, xl(n))
          ymin = min (ymin, yl(n))
          ymax = max (ymax, yl(n))
        enddo
c
c.... locate the grid points inside the block
c
c       hxy  = (1.66666666667d0)/float(100)
        hxy  = (660d0)/float(110)
        imin = int( (xmin - x0) / hxy + 1.9)
        imax = int( (xmax - x0) / hxy + 1.1)
        jmin = int( (ymin - y0) / hxy + 1.9)
        jmax = int( (ymax - y0) / hxy + 1.1)
c
        imin = max (imin, 1)
        imax = min (imax, maxi)
        jmin = max (jmin, 1)
        jmax = min (jmax, maxj)
c
        if (imin .gt. imax) goto 2000
        if (jmin .gt. jmax) goto 2000
c
        k    = 0
        do i = imin-1, imax-1
          do j = jmin-1, jmax-1
            k      = k + 1
            gxl(k) = x0 + i * hxy
            gyl(k) = y0 + j * hxy
        enddo
        enddo
c
        npgrd = (imax - imin + 1) * (jmax - jmin + 1)
        if (npgrd .gt. MAXPNT) goto 999
c
c.... project the data 
c
        call prj4q (xl,      yl,   gxl,     gyl,    
     &            soll,   npgrd, psoll,  lpsoll,     iel)
c
c.... copy the projected solution
c
        k = 0
        do i = imin, imax
          do j = jmin, jmax
            k = k + 1
            if (lpsoll(k) .eq. 1) then
              lpsol(i,j) = 1
              psol(i,j) = psoll(k)
            endif
          enddo
        enddo
c
c.... end of element loop
c
2000   continue
c
c.... return
c
      return
c
c.... error handling
c
999    call error ("projct  ", 
     &  "('Number of trial points exceeds MAXPNT = ',i6)", MAXPNT)
c
c.... end
c
      end
