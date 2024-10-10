        subroutine   velo_stiff_matrix_assemble
     &    (lm , idiag , stiff  , vlhsa , vlhsb, ivel)
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
      include 'common.h'
c
      dimension lm(numel,*), idiag(*), stiff(8,*), vlhsa(*), vlhsb(*) 
c
c ordering for skyline method
c
      index = 0
      do j = 1 , 8
        do i = 1 , 8              
          igadd = 0
          k = lm(ivel,j)
          if (k .ne. 0) then
            m = lm(ivel,i)
            if (m .ne. 0) then
              if(m .le. k  ) then 
                igadd  = idiag( k ) - k + m
              end if
            end if
          end if
          if (igadd .ne. 0) then
           vlhsa(igadd) = vlhsa(igadd) + stiff(i,j)
           vlhsb(igadd) = vlhsb(igadd) + stiff(j,i)
          end if
        enddo
      enddo
c
c.... return
c
       return
       end
