      subroutine untest(rhssym,rhsun,neq)
      implicit real*8 (a-h, o-z)
      dimension rhssym(*), rhsun(*)
 
      common /io    / iin,    igeom, iout , itsout , itout , imout, 
     &                irsin , irsout, icomp, igeoid
c     &     , ivisc , ivt  ,               itop  , ibot

      small = 1.0D-6
      write(iout,*) "begin untest"
      do 100 i=1, neq
      test = dabs(rhssym(i) - rhsun(i))
      if ( test .gt. small ) write(iout,*) i, rhssym(i), rhsun(i)
 100  continue
      write(iout,*) "end untest"
      return
      end
