      subroutine print(x,v,t,ntstep)
c
c.... program to print kinematic data
c
      use materials, only : tref  

      implicit double precision (a-h,o-z)
c
      include 'common.h'

      dimension x(2,*),v(2,*),t(*)

c.... temperature scale for temperature field for GMT

      write(itout,1000) ndof,nelx,nelz,numnp,ntstep,time
      write(itout,1500)
c

      do 200 n=1, numnp
        write(itout,2000) n,(x(i,n),i=1,ndof),(v(i,n),i=1,ndof),t(n),
     &                    (t(n)+tref(n)-T0/diff_T)
  200 continue


      return
c
 1000 format(5i10,f15.6)
 1500 format('___ node_____x1__________x2___________v1____________v2',
     &       '______tempature___temp+adiabat')
 2000 format(1x,i7,1x,2(1pe11.5,1x),4(1pe12.5,1x))
      end
