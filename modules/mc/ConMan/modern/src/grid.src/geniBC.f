      subroutine geniBC(id,ndof,numnp,neq,necho)
c
c.... program to read, generate and write boundary condition data
c        and establish equation numbers
c
c
      implicit double precision (a-h,o-z)
c
      dimension id(ndof,*)
c
        common /io    / iin,    igeom, iout , itsout , itout , imout, 
     &                  irsin , irsout, icomp, igeoid

      logical pflag
c
c
      call iclear(id,ndof*numnp)
      call igen(id,ndof)
c
      if (necho.eq.1) then
         nn=0
         do 200 n=1,numnp
         pflag = .false.
c
         do 100 i=1,ndof
         if (id(i,n).ne.0) pflag = .true.
  100    continue
c
         if (pflag) then      
            nn = nn + 1
            if (mod(nn,50).eq.1) then
              write(iout,1000) (i,i=1,ndof)
            endif
            write(iout,2000) n,(id(i,n),i=1,ndof)

         endif
  200    continue
      endif
c
c.... establish equation numbers
c
      neq = 0
c
      do 400 n=1,numnp
c
      do 300 i=1,ndof
      if (id(i,n).eq.0) then
         neq = neq + 1
         id(i,n) = neq
      else
         if (id(i,n) .eq. 1) then
           id(i,n) = 0
         else
           id(i,n) = -1
         endif
      endif
c
  300 continue
c
  400 continue
c
      return
c
 1000 format(' ',' n o d a l   b o u n d a r y   c o n d i t i o n   c o
     & d e s'///
     & 5x,' node no.',3x,6(6x,'dof',i1:)//)
 2000 format(6x,i7,5x,6(5x,i7))
c
      end
      
      subroutine igen(ia,m)
c
c.... program to read and generate integer nodal data
c
c        ia = input array
c         m = number of rows in ia
c         n = node number
c        ne = end node in generation sequence
c        ng = generation increment
c
c
      implicit double precision (a-h,o-z)
c
      dimension ia(m,*), ib(13)
      common /io    / iin,    igeom, iout , itsout , itout , imout, 
     &                irsin , irsout, icomp, igeoid
c     &           , ivisc , ivt  ,               itop  , ibot
c
  100 continue
1     read(iin,*,err=1,end=998) n,ne,ng,(ib(i),i=1,m)
      if (n.eq.0) return
      if (ng.eq.0) then
         ne = n
         ng = 1
      else
         ne = ne - mod(ne-n,ng)
      endif
c
      do 200 i=n,ne,ng
      call BCmove(ia(1,i),ib,m)
  200 continue
c
      go to 100
 998  call error('igen    ','end file',inn)
c
      end
c**** new *********************************************************************
      subroutine BCmove(ia,ib,n)
c
c.... program to move an integer array 
c
c
      implicit double precision (a-h,o-z)
c
      dimension ia(*),ib(*)
c
      do 100 i=1,n
      ia(i)=ib(i)
  100 continue
c
      return
      end
