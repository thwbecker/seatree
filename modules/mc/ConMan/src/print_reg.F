      subroutine print_reg(t)
c
c.... program to print kinematic data
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension t(111,101)
      common /io    / iin,igeom,iout ,itsout ,itout ,imout ,
     &                irsin ,irsout

      nn = 0
c
      do j=101,1,-1
        do i=1,111
          if(abs(t(i,j)) .lt. 1.0e-7) t(i,j) = 0.0d0
        enddo
      enddo
c
      do j=101,1,-1
          write(82,2000) (t(112-i,j),i=1,111)
c for GMT this works fine
c       do i=1,111
c         x = 660.0*(i-1)/110.
c         y = 600.0*(j-1)/100.
c         write(itout,*) x,y, t(i,j)
c       enddo
      enddo
c
c
      return
c
 2000 format(111(1pe14.7,1x))
      end
