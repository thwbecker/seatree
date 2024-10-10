      subroutine genBC(f,ndof,numnp,j,necho, x, idv)
c
c.... program to read, generate and write nodal input data
c
c        f(ndof,numnp) = prescribed forces/kinematic data 
c                             =       l velocities (j=1)
c                             =        temperatures(j=2)
c
      use sub_bench
      implicit double precision (a-h,o-z)
c
      logical lzero
      dimension f(ndof,numnp), x(2,numnp), idv(ndof,numnp)

      common /const /zero  , pt25,    pt33,   pt5,     pt66,   one,   
     &               onept5, two ,   three,  four,   sixten,  eps7
      common /io    / iin,    igeom, iout , itsout , itout , imout, 
     &             irsin , irsout, icomp, igeoid



      
      call clear(f,numnp*ndof)
c
      call genfl(f,ndof)
      call ztest(f,ndof*numnp,lzero)
c
c added to use Batchelor Boundary Conditions
c
      if (ndof .gt. 1) then
      do n=1,numnp
        do i=1,ndof
          if (idv(i,n) .eq. -1) then
             call batchelor(x(1,n),x(2,n),f(1,n),f(2,n),
     &            sb_box_height, sb_plate_thickness, sb_uplate)
c             write(iout,*) i,n,idv(i,n), idv(2,n)
              idv(i,n) = 0
          endif
        enddo
      enddo
      endif
c
      if (necho.eq.1) then
c
         if (lzero) then
c
            if (j.eq.1) then
               write(iout,2000)
c
              endif
            if (j.eq.2) then
               write(iout,3000)
c
              endif
         else
c
             if (j.eq.1) 
     &    call printd(' B o u n d a r y    V e l o c i t i e s  ',
     &                  f,ndof,numnp,0,ZERO)
c
            if (j.eq.2) 
     &    call printd(' B o u n d a r y  T e m p e r a t u r e  ',
     &                  f,ndof,numnp,0,ZERO)
c
         endif
      endif
c
c
      return
 2000 format(' '//,' there are no nonzero Boundary velocities')
 3000 format(' '//,' there are no nonzero Boundary temperature   ')
      end

      subroutine ztest(a,n,lzero)
c
c.... program to determine if an array contains only zero entries
c
      implicit double precision (a-h,o-z)
c
      dimension a(*)
      common /const / zero  ,   pt25,     pt33,   pt5,     pt66,   one,
     &               onept5 ,   two ,    three,  four,   sixten,  eps7
      logical lzero
c
      lzero = .true.
c
      do 100 i=1,n
      if (a(i).gt.eps7) then
         lzero = .false.
         return
      endif
  100 continue
c
      return
      end 
