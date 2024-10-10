      subroutine printd(name,dva,ndof,numnp,ntstep,time)
c
c.... program to print kinematic data
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      logical lzero,lskip
      dimension name(11),dva(ndof,1)
      common /io    / iin, igeom, iout , itsout , itout , imout , 
     &                irsin , irsout
c
      nn = 0
      lskip = .true.
c
      do 100 n=1,numnp
      call ztest(dva(1,n),ndof,lzero)
      if (.not.lzero) then
         nn = nn + 1
         if (mod(nn,50).eq.1) then
           write(iout,1000) name,ntstep,time,(i,i=1,ndof)
      endif
         write(iout,2000) n,(dva(i,n),i=1,ndof)
         lskip = .false.
      endif
  100 continue
c
      if (lskip) then
         write(iout,1000) name,ntstep,time,(i,i=1,ndof)
         write(iout,3000)  

      endif
c
      return
c
 1000 format(' ',11a4//5x,
     &' step number = ',i10//5x,
     &' time        = ',1pe10.3///5x,
     &' node no.',6(13x,'dof',i1,:)/)
 2000 format(6x,i5,10x,6(1pe15.8,2x))
 3000 format(' ',//,' there are no nonzero components')
      end
