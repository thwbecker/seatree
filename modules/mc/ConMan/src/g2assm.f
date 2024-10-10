        subroutine   g2assm(lm , idiag , stiff  , vlhs , ivel )
c----------------------------------------------------------------------
c
c   This subroutine assembles the local LHS mass matrices into the 
c global upper  matrix.
c
c input:
c    lm  (numel,8)         : element mapping array
c    idiag  (nEGdf)        : index of diag. entries of LHS Matrix
c    stiff (36)            : element stiff matrix 
c    nenl                  : number of local nodes for this element block
c
c output:
c    vlhs (nsize)          : LHS matrix - upper triangle
c
c
c----------------------------------------------------------------------
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
        include 'common.h'
c
        dimension lm(numel,*), idiag(*), stiff(*), vlhs(*) 
c
c ordering for skyline method
c
      if (isky .eq. 1) then
        index = 0
        do 500 j = 1 , 8
        do 400 i = 1 , j
        index = index + 1              
        igadd = 0
        k = lm(ivel,j)
        if (k .ne. 0) then
          m = lm(ivel,i)
          if (m .ne. 0) then
            if(k .ge. m  ) then 
              igadd  = idiag( k ) - k + m
            else
              igadd  = idiag( m ) - m + k
            end if
          end if
        end if
       if (igadd .ne. 0) then
          vlhs(igadd)   = vlhs(igadd) + stiff(index)
        end if
300     continue
400     continue
500     continue
      elseif (isky .eq. 0) then
c
c.... ordering for vector factor/solve
c
        index = 0
        do 900 j = 1 , 8
        do 800 i = 1 , j
        index = index + 1              
        igadd = 0
        k = lm(ivel,j)
        if (k .ne. 0) then
          m = lm(ivel,i)
          if (m .ne. 0) then
            if(k .le. m  ) then 
              igadd = idiag( k ) - k + m
            else
              igadd = idiag( m ) - m + k
            end if
          end if
        end if
        if (igadd .ne. 0) then
          vlhs(igadd) = vlhs(igadd) + stiff(index)
        end if
700     continue
800     continue
900     continue
      end if
c
c.... return
c
       return
       end
