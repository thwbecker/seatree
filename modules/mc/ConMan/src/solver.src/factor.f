      subroutine factor(a,idiag,neq,temp)
c
c.... program to perform Crout factorization: a = u(transpose) * d * u
c
c        a(i):  coefficient matrix stored in compacted column form;
c               after factorization contains d and u
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      dimension a(*),idiag(*),temp(*)
      common /const /zero  , pt25,    pt33,   pt5,     pt66,   one,
     &               onept5, two ,   three,  four,   sixten,  eps7
      jj = 0 
c
      do 300 j=1,neq
c
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
c
      if (jcolht.gt.2) then
c  
c....... for column j and i.le.j-1, replace a(i,j) with d(i,i)*u(i,j)
c
         istart = j - jcolht + 2
         jm1    = j - 1
         ij     = jjlast + 2
         ii     = idiag(istart-1)
c
         do 100 i=istart,jm1
c  
         iilast = ii
         ii     = idiag(i) 
         icolht = ii - iilast
         jlngth = i - istart + 1  
         length = min0(icolht-1,jlngth)  
         if (length.gt.0)  
     &      a(ij) = a(ij) - coldot(a(ii-length),a(ij-length),length) 
         ij = ij + 1   
  100    continue 
c   
      endif
c 
      if (jcolht.ge.2) then
c 
c....... for column j and i.le.j-1, replace a(i,j) with u(i,j);
c           replace a(j,j) with d(j,j).
c
         jtemp = j - jj
c 
         SUM=zero
c$dir no_recurrence
         do 200 ij=jjlast+1,jj-1
            temp(ij-jjlast)  = a(ij) 
            a(ij) = a(ij)/a( idiag(jtemp + ij) )
            SUM = SUM + temp(ij-jjlast)*a(ij)
  200    continue 
         a(jj) = a(jj) - SUM
c
      endif
c
  300 continue
c
      return
      end
