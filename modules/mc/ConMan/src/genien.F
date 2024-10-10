#ifdef IMPLICIT
         subroutine genien (ien   , idiag ,
     &                      mat   , lmv    , lmt   , id    ,
     &                      idiagt, idt)
#else
         subroutine genien (ien   , idiag ,
     &                      mat   , lmv    , lmt   , id    )
#endif
c
c----------------------------------------------------------------------
c this program reads and generates the ien,idaig arrays
c
c output:
c    ien    (numel,nen)   : input from GEOM.DAT as element connectivity 
c    idiag  (neqv)        : index of Diag. entry of Macro Elem. LHS 
c    lmv    (numel,8)     :  mapping array for velocity
c    lmt    (numel,4)     :  mapping array for temperature
c    nsize                : Total size of Macro Elem. LHS Matrix
c
c
c----------------------------------------------------------------------
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
#include<common.h>
c
      common /temp1 / itemp(8)
c
      dimension ien (numel,*) ,
     &          idiag(*)      , mat(*)   , 
     &          lmv(numel,*)  , lmt(numel,*)  , id(2,*)  
#ifdef IMPLICIT
      dimension idiagt(*)     , idt(*)
#endif
c.... -------------------->  Generate IEN  <--------------------
c
c.... read element nodal data card
c
100      continue
      call iclear (itemp, nen)
      n = 0
1     read (igeom,*,err=1,end=998) n, ng, mat_no ,(itemp(i),i=1,nen)
      if (n. eq. 0) goto 500
      mat(n) = mat_no
      do 200 i = 1, nen
        ien(n,i) = itemp(i)
200   continue
      if (ng .eq. 0) goto 100
c
c.... initialize the generation data and input the changes
c
      n1 = 1
      n2 = 1
      n3 = 1
      ince1 = 1
      ince2 = 0
      ince3 = 0
      inc1 = itemp(2) - itemp(1)
      inc2 = itemp(3) - itemp(2)
      inc3 = itemp(5) - itemp(1)
c
2      read (igeom,*,err=2,end=998) n1, ince1, inc1, n2, ince2, inc2
c     &           inc2, n3, ince3, inc3
c
      if (ince2 .eq. 0) ince2 = n1
      if (ince3 .eq. 0) ince3 = n1 * n2
c
c.... generate ien
c
      do 400 i=0,n1-1
      do 400 j=0,n2-1
      do 400 k=0,n3-1
          iel  = n + i * ince1 + j * ince2 + k * ince3
          if (iel .gt. numel) call error ('genien  ','elem num',iel)
          incr = i * inc1 + j * inc2 + k * inc3
          mat(iel) = mat_no 
          do 300 nn=1,nen
            ien(iel,nn) = itemp(nn) + incr
300       continue
400   continue
      goto 100
c
c.... echo the results
c
500   continue
      if (necho .eq. 1) then
        do 600 iel = 1, numel
          if (mod(iel,50).eq.1) write (iout,5000) ititle,(i,i=1,nen)
          write (iout,5101) iel, mat(iel), (ien(iel,i),i=1,nen)
600     continue
      endif
c
c.... create the lm arrays from array
c
      do 700 ivel  = 1, numel
        lmv(ivel,1)  = id(1, ien(ivel,1) )
        lmv(ivel,2)  = id(2, ien(ivel,1) )
        lmv(ivel,3)  = id(1, ien(ivel,2) ) 
        lmv(ivel,4)  = id(2, ien(ivel,2) )
        lmv(ivel,5)  = id(1, ien(ivel,3) )
        lmv(ivel,6)  = id(2, ien(ivel,3) )
        lmv(ivel,7)  = id(1, ien(ivel,4) )
        lmv(ivel,8)  = id(2, ien(ivel,4) )
#ifdef IMPLICIT
        lmt(ivel,1)  = idt(  ien(ivel,1) )
        lmt(ivel,2)  = idt(  ien(ivel,2) )
        lmt(ivel,3)  = idt(  ien(ivel,3) )
        lmt(ivel,4)  = idt(  ien(ivel,4) )
#else
        lmt(ivel,1)  = ien(ivel,1)
        lmt(ivel,2)  = ien(ivel,2)
        lmt(ivel,3)  = ien(ivel,3)
        lmt(ivel,4)  = ien(ivel,4)
#endif
700   continue
c
c.... if wrap around condition let the left side nodes of
c     the first elements equal the right side nodes of the last
c     elements
c
      if (nwrap .ne. 0) then
         do 800 iel = 1, nwrap
           iwrap = nel - nwrap + iel
           lmv(iwrap,3) = lmv(iel,1)          
           lmv(iwrap,4) = lmv(iel,2)          
           lmv(iwrap,5) = lmv(iel,7)          
           lmv(iwrap,6) = lmv(iel,8)          
           lmt(iwrap,2) = lmt(iel,1)          
           lmt(iwrap,3) = lmt(iel,4)          
800      continue
      end if
c
c     find column hts
c     
      if (isky .eq. 1) then
        do 2600 k=1,numel
        min = neqv 
c 
        do 2100 j=1,8
c 
          num = lmv(k,j)
          if (num.gt.0) min = min0(min,num)
 2100   continue 
c
        do 2400 j=1,8
c
          num = lmv(k,j)
          if (num.gt.0) then
             m = num - min
             if (m.gt.idiag(num)) idiag(num) = m   
          endif 
c  
 2400   continue
c  
 2600   continue 
c 
c
c.... convert column height to diagonal index
c
        n = 1
        idiag(1) = 1 
        if (neqv .eq. 1) return
c
        do 2700 i=2,neqv
        idiag(i) = idiag(i) + idiag(i-1) + 1
 2700   continue
c
      nsize = idiag(neqv)
      else
c
        do 3600 k=1,numel 
        max = 0   
c   
        do 3100 j=1,8
c                 
        num = lmv(k,j)
        if (num.gt.0) max = max0(max,num)
 3100 continue 
c            
        do 3400 j=1,8 
c    
        num = lmv(k,j)
        if (num.gt.0) then
         m = max - num  
         if (m.gt.idiag(num)) idiag(num) = m
        endif  
c
 3400 continue
c    
 3600 continue
c
c
c.... convert column height to diagonal index
c
        long = 0
        do 3700 i = 1, neqv
          long = max0(long,idiag(i))
 3700   continue
        idiag(1) = 1
        if (neqv .eq. 1) return
        do 3800 i = 2, neqv
          ltemp = i + long - 1
          if (ltemp .gt. neqv ) then
            ntemp = neqv - i + 1
          else
            ntemp = long
          end if
          idiag(i) = ntemp + idiag(i-1) + 1
 3800   continue
        nsize= idiag(neqv)
      end if
#ifdef IMPLICIT
c
c find idiag for temperature equation
c
      do 5600 k=1,numel
        min = neqt 
c 
        do 5100 j=1,4
c 
          num = lmt(k,j)
          if (num.gt.0) min = min0(min,num)
 5100   continue 
c
        do 5400 j=1,4
c
          num = lmt(k,j)
          if (num.gt.0) then
             m = num - min
             if (m.gt.idiagt(num)) idiagt(num) = m   
          endif 
c  
 5400   continue
c  
 5600   continue 
c 
c
c.... convert column height to diagonal index
c
        n = 1
        idiagt(1) = 1 
        if (neqt .eq. 1) return
c
        do 5700 i=2,neqt
        idiagt(i) = idiagt(i) + idiagt(i-1) + 1
 5700   continue
c
      nsizet = idiagt(neqt)
#endif
c
c
c.... return
c
      return
c
c.... end of file error handling
c
998      call error ('genien  ','end file',igeom)
c
5000      format(' ',80a1,//,
     &  ' E l e m e n t   N o d a l   C o n n e c t i v i t y',//,
     &  '   Elem ',/,
     &  '  Number  ',5x,' mat ',8x,8('Node',i1,:,3x))
5101      format(2x,i5,6x,i5,6x,8i8)
      end
