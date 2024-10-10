c     box height (600)
c     plate thickness (50)

      subroutine fault(v, xl, uq, vq, is_a_fault,box_height,
     &     plate_thickness)
c
c Written specifically for the geometry of the subduction
c benchmark.
c
c     implicit double precision (a-h,o-z) 
      dimension v(2,*),xl(2,4), uq(5), vq(5)
c
#include<common.h>
      double precision poff,box_height,plate_thickness

      poff = box_height - plate_thickness
c
c  if we are in the top 50 km / pt
c
      if (xl(2,1) .gt. poff) then
c
c  if we are in an element along the diagnoal (x = z) which contains
c  the fault
c
        if(abs(xl(1,2)-(box_height - xl(2,2))) .lt. 0.01) then
c
c Gauss point 1 is below the fault, gauss points 2, 4 are on the fault
c and gauss point 3 is above the fault, v(*,nelz+1) should be the top node 
c at (0,600) and is Vslab in the geometery used here. This is very specific.
c Rather than interpolate, we are setting the gauss points to be either Vslab
c or zero.
c
          is_a_fault = 1
          uq(1) = v(1,nelz+1)
          vq(1) = v(2,nelz+1)
          uq(2) = v(1,nelz+1)
          vq(2) = v(2,nelz+1)
          uq(3) = 0.0d0
          vq(3) = 0.0d0 
          uq(4) = v(1,nelz+1)
          vq(4) = v(2,nelz+1)
        endif
      endif
c
c that's it
c
      return
      end
