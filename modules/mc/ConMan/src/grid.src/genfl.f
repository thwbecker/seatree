      subroutine genfl (dest, ndest)
c
c----------------------------------------------------------------------
c
c  This routine reads and/or generates floating-point nodal data.
c
c input:
c    ndest                 : number of rows in dest (le.6)
c
c output:
c    dest (ndest,numnp)    : the destination array
c
c----------------------------------------------------------------------
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      include '../common.h'
c
      dimension temp(12,20), ninc(3), inc(3)
      dimension dest(ndest,1)
c
c.... read the first card (if node.eq.0 terminate the generation)
c
100   continue
      node = 0
      call clear (temp(1,1),12*20)
1     read (igeom,*,err=1,end=998) node,numgp,(temp(i,1),i=1,ndest)
c
      if (node .eq. 0) goto 300
c
      call smove (dest(1,node), temp, ndest)
      if (numgp .eq. 0) go to 100
c
c.... read the generation cards
c
      do 200 j=2,numgp
        call clear (temp(1,j),ndest)
2       read (igeom,*,err=2,end=998) m, mgen, (temp(i,j),i=1,ndest)
        if (mgen .ne. 0) call move (dest(1,m), temp(1,j), ndest)
200   continue
c
c.... set up the increments and generate the data
c
      call iclear (ninc,3)
      call iclear (inc,3)
3     read (igeom,*,err=3,end=998) (ninc(i),inc(i),i=1,ndest)
      call genfl1 (temp, ninc, inc, node, numgp, ndest, dest)
c
c... return to the beginning of the loop for more data
c
      go to 100
c
c.... return
c
300   return
c
c.... end of file error handling
c
998   call error ('genfl   ','end file',igeom)
c
      end
