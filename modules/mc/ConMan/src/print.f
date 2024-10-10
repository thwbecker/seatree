      subroutine print(x,v,t,ndof,nelx,nelz,numnp,ntstep,time)
c
c.... program to print kinematic data
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension x(ndof,1),v(ndof,1),t(1)
      common /io    / iin,igeom,iout ,itsout ,itout ,imout ,
     &                irsin ,irsout

      nn = 0
c
      write(itout,1000) ndof,nelx,nelz,numnp,ntstep,time
      write(itout,1500)
      do 100 n=1,numnp
      write(itout,2000) n,(x(i,n),i=1,ndof),(v(i,n),i=1,ndof),t(n)
  100 continue
c
c
      return
c
 1000 format(5i10,f10.6)
 1500 format('  node   x1        x2           v1             v2',
     &       '         tempature    ')
 2000 format(1x,i7,1x,2(1pe11.5,1x),5(1pe12.5,1x))
      end
