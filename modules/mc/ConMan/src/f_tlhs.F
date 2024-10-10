      subroutine f_tlhs( tlhsa,   tlhsb,     shl,
     &                    shdx,   shdy,    
     &                   lmt,     trhs,    idiag,   ivel,
     &                   det)     
c---------------------------------------------------------------------
c
c  This routine calculates the L.H.S. temperature matrix
c  
c----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
c
#include<common.h>
#include<materials.h>
c
c
      dimension  shdx(4,*), shdy(4,*), shl(4,*)     ,
     &           det(*) , idiag(*)       , lmt(numel,*)   ,
     &           tlhsa(*)     , tlhsb(*)       , trhs(*)        
c
      common /temp1 / el_rhs(4),blkmhu , pg(4,5) ,
     &                vl(8)    , uq(5) , vq(5)   ,
     &                conv     ,tmass  , adiff   ,
     &                unorm    ,eta    , xse     ,
     &                blkdif   ,uxse   , ueta    ,
     &                tlhsl(4,5),tl(8) , 
     &                m      ,k      , igadd  
c
      do i=1,4
        do j=1,4
          tlhsl(i,j) = zero
        enddo
      enddo
c
      do intp = 1 , 4
c
c add advection part of l.h.s.
c   template:
c        tlhsl(i,j)  = tlhsl(i,j) + pg(i,l)*
c    &              (shdx(j,l)*uq(l) + shdy(j,l)*vq(l)) 
c    &              * det(l)
c
        do nodej = 1, 4
          conv =  (shdx(nodej,intp)*uq(intp) 
     &           + shdy(nodej,intp)*vq(intp)) * det(intp) 
          do nodei = 1, 4
             tlhsl(nodei,nodej)  = tlhsl(nodei,nodej) + 
     &                 pg(nodei,intp)* conv
          enddo
        enddo
c
c add diffusion part of l.h.s.
cc  template:
cc       tlhsl(i,j)  = tlhsl(i,j) + (shdx(ivel,i,l)*shdx(ivel,j,l) 
cc   &              + shdy(ivel,i,l)*shdy(ivel,j,l)) * det(ivel,l)*blkdif
c
        blkmhu = dmhu(mat(ivel))
        do nodej = 1, 4
          do nodei = 1, 4
            tlhsl(nodei,nodej) = tlhsl(nodei,nodej)
     &        +(shdx(nodej,intp)*shdx(nodei,intp)
     &        + shdy(nodej,intp)*shdy(nodei,intp))
     &        *det(intp)*blkdif
          enddo
          el_rhs(nodej) = el_rhs(nodej)
     &            +shl(nodej,intp)*blkmhu*det(intp)
        enddo
c
c end integration point loop
c
      enddo
c
c update temperature r.h.s. residual
c
      do node = 1, 4
        el_rhs(node)  =  el_rhs(node)
     &       - tlhsl(node,1)*tl(1) - tlhsl(node,2)*tl(2)
     &       - tlhsl(node,3)*tl(3) - tlhsl(node,4)*tl(4)
      enddo
c
c add mass term
c 
#ifdef PICARD
#else
      do j = 1 , 4
        do i = 1 , 4
          tlhsl(i,j)  = tlhsl(i,j)*alpha*dt
        enddo
      enddo
c
      do intp = 1,4
        do nodej = 1,4
          do nodei = 1,4
          tmass = (shl(nodej,intp) )*det(intp)
          tlhsl(nodei,nodej) = tlhsl(nodei,nodej) + pg(nodei,intp)*tmass
          enddo
        enddo
      enddo
c
#endif
c
c assemble nonsymmetric ta and tc
c ordering for skyline method
c
        do j = 1 , 4
          do i = 1 , 4
            igadd = 0
            k = lmt(ivel,j)
            if (k .ne. 0) then
              m = lmt(ivel,i)
              if (m .ne. 0) then
                if(m .le. k ) then 
                  igadd  = idiag( k ) - k + m
                end if
              end if
            end if
            if (igadd .ne. 0) then
            tlhsa(igadd) = tlhsa(igadd) + tlhsl(i,j)
            tlhsb(igadd) = tlhsb(igadd) + tlhsl(j,i)
          end if
        enddo
      enddo
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      return
      end  
