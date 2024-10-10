      subroutine vfactor(a,idiag,neq,temp,nsizea)
c
c.... program to perform Crout factorization: a = u(transpose) * d * u
c
c        a(i):  coefficient matrix stored in compacted column form;
c               after factorization contains d and u
c
c
      implicit double precision (a-h,o-z)
c
      dimension a(nsizea),idiag(neq),temp(neq)
      common /const / zero ,   pt25,     pt33,   pt5,    pt66,   one,
     &   onept5,   two ,    three,  four,  sixten,  eps7
c
c
      ndiag = idiag(neq)
      do 1000 k = 1 , neq - 1
          kdiag = idiag(k)
          akdiag = one / a(kdiag)
          k1diag = idiag(k+1)
          klenth = k1diag - kdiag - 1
          koff = 0
          if( klenth .ge. 1 ) then
c
c$dir no_recurrence
c
!$omp parallel do if (klenth .ge. 1000)
              do 100 i = 1 , klenth
                 temp(i) = a(kdiag + i ) * akdiag
100           continue
!$omp end parallel do 
!$omp parallel do if (neq - klenth + 1 .ge. 1000)
              do 101 i = klenth + 1, neq
                 temp(i) = zero 
101           continue
!$omp end parallel do 
c
              do 300 j = 1 ,klenth 
                if( k + j .lt. neq ) then
                   const = a(kdiag + j )
                   if ( const .ne. zero ) then
                     jdiag = idiag(j + k ) 
                     j1diag = idiag(j + k + 1)
                     jlenth = j1diag - jdiag 
                     if ( klenth .le. jlenth - 1) then
                       lenth = jlenth  - koff - 1 
                     else
                       lenth = jlenth
                     end if
c$dir no_recurrence
!$omp parallel do if (lenth .ge. 1000)
                     do 200 i = 1 , lenth
                        a(jdiag + i  - 1) = a(jdiag + i  - 1) -
     &                                     const * temp(i + koff)
200                  continue
!$omp end parallel do 
                   end if
                   koff = koff + 1
c
                else
                 a(ndiag) = a(ndiag) - temp(klenth) * a(kdiag+klenth)
                end if
300           continue
c    
c$dir no_recurrence
c
!$omp parallel do if (klenth .ge. 1000)
              do 400 i = 1 , klenth
                  a(kdiag + i ) = temp(i) 
400           continue
!$omp end parallel do 
          end if
c
1000   continue
c
       return
       end
c**** new **********************************************************************
      subroutine vback(a,b,idiag,neq,nsizea)
c
c.... program to perform forward reduction and back substitution
c
c
      implicit double precision (a-h,o-z)
c
      dimension a(nsizea),b(neq),idiag(neq)
      common /const / zero ,   pt25,     pt33,   pt5,    pt66,   one,
     &               onept5,   two ,    three,  four,  sixten,  eps7
c
c.... forward reduction
c
      do 200 k = 1 , neq - 1
        const = b(k)
        if ( const .ne. zero ) then
        k1diag = idiag(k+1)
        kdiag = idiag(k)
        klenth = k1diag - kdiag - 1
c
c$dir   no_recurrence
!$omp parallel do if (klenth .ge. 1000)
           do 100 m = 1 , klenth
              b(m + k) = b(m + k ) - const * a( kdiag + m )
100        continue
!$omp end parallel do 
c 
      end if
200   continue
c
c.... diagonal scaling
c
c$dir   no_recurrence
!$omp parallel do if (neq .ge. 1000)
      do 350 k = 1 , neq
         b(k) = b(k) / a(idiag(k))
350   continue
!$omp end parallel do 
c
c.... backsubtitution
c
      do 400 k = neq - 1 , 1 , -1
        kdiag = idiag(k)
        k1diag = idiag(k+1)
        lenth = k1diag - kdiag -1
        ndiag = kdiag
c$dir   no_recurrence
         psum = b(k)
!$omp parallel do if (lenth .ge. 1000) reduction (+: psum)
             do 300 j = 1 , lenth
                psum =  psum -  a(kdiag + j ) * b(k + j)
300          continue
!$omp end parallel do 
         b(k) = psum
400   continue
      return
      end
