      subroutine prtstr_reg(press, iflag)
c
c.... program to print stress data
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension press(111,101)
c
      common /io    / iin,igeom,iout ,itsout ,itout ,imout ,
     &                irsin ,irsout

      nn = 0
c
c     pscale = press(2,2)
      pscale = press(110,2)
      if (iflag .eq. 1) then
      do j=101,1,-1
        do i=1,111
c         press(i,j) = press(i,j) - pscale
        enddo
      enddo
      endif
 
      do j=101,1,-1
c         write(imout,2000) (press(112-i,j),i=1,111)
c for GMT this works fine
        do i=1,111
          x = 660.0-660.0*(i-1)/110.
          y = 600-600.0*(j-1)/100.
          write(imout,*) x,y, press(i,j)
        enddo
      enddo
c
c
      return
c
 2000 format(111(1pe14.7,1x))
      end
