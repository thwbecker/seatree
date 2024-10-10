      subroutine EGlib (task)
c
c----------------------------------------------------------------------
c
c  This routine is the Element Group LIBrary routine. This is the only 
c link between the global driver and the element group routines.
c
c input:
c  mnpar      : Pointer to the parameters of the element group
c  task       : task to be performed
c
c----------------------------------------------------------------------
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      include 'common.h'
c
      character*8 task
c
c.... find the element type
c
      goto (200, 300) ntype-1
c
      call error ('EGlib   ', 'ntype   ', ntype)
      return
c
c.... 2D Element Group
c
200   call EG2 (task)
      return
c
c.... 3D Element Group
c
300   continue
      return
c
      end
