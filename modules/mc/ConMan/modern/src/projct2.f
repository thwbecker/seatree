      subroutine projct (sol, psol, ipsol ) 
c
c-----------------------------------------------------------------------
c
c  This routine projects the solution from an unstructured mesh to
c a structured grid.
c
c input:
c  sol        	: solution on unstructured mesh
c  npsol        : array size of grid solution
c
c output:
c  psol        	: projected solution
c  lpsol        : flag to indicate a projected solution
c
c Farzin Shakib, Spring 1989.
c-----------------------------------------------------------------------
c
       use fe
       use variables

       implicit double precision (a-h,o-z)
       parameter (ipsolx = 301, ipsolz = 101)
       parameter (ipmax  = ipsolx * ipsolz)

       include 'common.h'
c
       dimension    sol(*), psol(ipsolx,ipsolz),lpsol(ipsolx,ipsolz),
     &              ipsol(ipsolx,ipsolz)   
       
c
       dimension    xl(4),  yl(4),  gxl(ipmax), gyl(ipmax),
     &            soll(4),  psoll(ipmax), lpsoll(ipmax)
c
       integer maxi, maxj

c.... initialization
        do i = 1, ipsolx
          do j = 1, ipsolz
            lpsol(i, j) = 0
           enddo
        enddo
        do i = 1, ipmax
          gxl(i) = zero ; gyl(i) = zero ; psoll(i) = zero
          lpsoll(i) = 0
        enddo
        do i = 1, 4
          xl(i) = zero ; yl(i) = zero ; soll(i) = zero
        enddo

c
c.... setup some constants
c
      x0   = x(1,1)
      y0   = x(2,1)
c
      maxi = ipsolx
      maxj = ipsolz            

c      maxi = 65
c      maxj = 65            


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
        do n = 1, 4
          node     = ien(iel,n)
          xl(n)    = x(1,node)
          yl(n)    = x(2,node)
          soll(n)  = sol(node)
        enddo
c
c.... compute the rectangular block surrounding the element
c
        xmin = xl(1)
        xmax = xl(1)
        ymin = yl(1)
        ymax = yl(1)
c
        do n = 2, nen
          xmin = min (xmin, xl(n))
          xmax = max (xmax, xl(n))
          ymin = min (ymin, yl(n))
          ymax = max (ymax, yl(n))
        enddo
c
c.... locate the grid points inside the block
c
c       hxy  = (1.66666666667d0)/float(100)

c.... for lectangular box (length of bottom side is 660)

c        hxy  = x(1,numnp)/float(nelx)
c        imin = int( (xmin - x0) / hxy + 1.9)
c        imax = int( (xmax - x0) / hxy + 1.1)
c        jmin = int( (ymin - y0) / hxy + 1.9)
c        jmax = int( (ymax - y0) / hxy + 1.1)

c.... for square box (length of side is 1.0)

c        hxy  = xsize/float(nelx)
        hxy  = x(1,numnp)/float((ipsolx-1))

        imin = int( (xmin - x0) / hxy + 1.1)
        imax = int( (xmax - x0) / hxy + 1.1)
        jmin = int( (ymin - y0) / hxy + 1.1)
        jmax = int( (ymax - y0) / hxy + 1.1)
        
c
c        write(*,*) imin, imax, jmin, jmax
        
        imin = max (imin, 1)
        imax = min (imax, maxi)
        jmin = max (jmin, 1)
        jmax = min (jmax, maxj)
c
        if (imin .gt. imax) goto 2000
        if (jmin .gt. jmax) goto 2000
c
c        write(*,*) 'did'
        
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
        if (npgrd .gt. ipsolx) goto 999
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
              ipsol(i,j) = int(255*psoll(k))
c              write(*,*) psol(64,j), ipsol(64,j)
            endif
          enddo
        enddo
        
c          if ((time .gt. 0.009000) .and. (time .lt. 0.0090020)) then
c           do j = 1, 65
c             write(iout,2001) (psol(i,j), i = 1,601,1)
c           enddo
c           write(iout, 2001)
c          endif
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
     &  "('Number of trial points exceeds ipsolx = ',i6)", ipsolx)
c2001  format(601(f15.5))
c
c.... end
c
      end
