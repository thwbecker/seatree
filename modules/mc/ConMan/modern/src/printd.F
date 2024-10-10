      subroutine printd(name,dva,ndof,numnp,ntstep,time)
c
c.... program to print kinematic data
c
      implicit double precision (a-h,o-z)
c
      logical lzero,lskip
      character*1 name(41)
      dimension dva(ndof,numnp)
      common /io    / iin,    igeom, iout , itsout , itout , imout, 
     &                  irsin , irsout, icomp, igeoid
c     &                , ivisc , ivt  ,  ivin  , ivout
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
 1000 format(' ',41a1//5x,
     &' step number = ',i10//5x,
     &' time        = ',1pe10.3///5x,
     &' node no.',6(13x,'dof',i1,:)/)
 2000 format(6x,i5,10x,6(1pe15.8,2x))
 3000 format(' ',//,' there are no nonzero components')
      end
