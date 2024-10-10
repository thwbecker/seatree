      subroutine untest(rhssym,rhsun,neq)
      implicit real*8 (a-h, o-z)
      dimension rhssym(1), rhsun(1)
      common /io    / iin,igeom,iout ,itsout ,itout ,imout ,
     &                irsin ,irsout
      small = 1.0e-6
      write(iout,*) "begin untest"
      do 100 i=1, neq
      test = abs(rhssym(i) - rhsun(i))
      if ( test .gt. small ) write(iout,*) i, rhssym(i), rhsun(i)
 100  continue
      write(iout,*) "end untest"
      return
      end
