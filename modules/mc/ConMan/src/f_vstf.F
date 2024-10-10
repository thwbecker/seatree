      subroutine f_vStf(x     , v     , t     , idiag , 
     &                  vlhs  , vbcr  , vrhs  )
c
c----------------------------------------------------------------------
c
c  This routine calculates the L.H.S. stiffness matrix 
c  for the Implicit Element Group with Macro-Element
c  preconditioner.  The LHS  matrix are statically condensed 
c  with respect to the internal nodes. 
c  **** note  if there is only one group then all nodes are internal. 
c             
c----------------------------------------------------------------------
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
#include<common.h>
#include<materials.h>
#include<fe.h>
#include<bc.h>
c
      dimension  x(nsd,*)     , v(ndof,*)   , t(*)     ,
     &           idiag(*)     , vlhs(*)     ,
     &           vbcr(*)      , vrhs(*)     ,
     &           shdx(4,5)    , shdy(4,5)   , det(5)
c
c
      logical lfirst
      common /temp1 / stiff(36),blklam,blkvis,power,vbcl(8),tq(5),
     &                tl(8)    , evisc(5), xl(2,4)
c
      if (isky .ne. 2) then
        do i = 1 , nsize
         vlhs(i) = zero
        enddo
      endif
      do i=1, neqv
         vrhs(i) = zero
         vbcr(i) = zero
      enddo
c
c.... loop over the element blocks
c
      do 1000 ivel = 1 , numel
c
c.... set up the parameters
c
      call rheol (x , v , t , ivel , evisc)
c
      do n = 1, 4
        vbcl(2*n-1) = vbc(1, ien(ivel,n) )
        vbcl(2*n)   = vbc(2, ien(ivel,n) )
        tl(n)       = t(  ien(ivel,n) ) 
        xl(1,n)     = x(1, ien(ivel,n) ) 
        xl(2,n)     = x(2, ien(ivel,n) ) 
      enddo
c
c.... call the global shape function
c
        call genshg (shdx, shdy , det , xl , ivel)
c
c
c... form the element stiffness for this block
c
c
c....form the lamda stiffness 
c
        blklam = alam(mat(ivel))*evisc(5)*det(5)*four
c       blklam = alam(mat(ivel))*det(5)*four
c
        stiff(1) = shdx(1,5) * shdx(1,5) * blklam
        stiff(2) = shdx(1,5) * shdy(1,5) * blklam
        stiff(3) = shdy(1,5) * shdy(1,5) * blklam
        stiff(4) = shdx(1,5) * shdx(2,5) * blklam
        stiff(5) = shdx(2,5) * shdy(1,5) * blklam
        stiff(6) = shdx(2,5) * shdx(2,5) * blklam
        stiff(7) = shdx(1,5) * shdy(2,5) * blklam
        stiff(8) = shdy(1,5) * shdy(2,5) * blklam
        stiff(9) = shdx(2,5) * shdy(2,5) * blklam
        stiff(10) = shdy(2,5) * shdy(2,5) * blklam
        stiff(11) = shdx(3,5) * shdx(1,5) * blklam
        stiff(12) = shdx(3,5) * shdy(1,5) * blklam
        stiff(13) = shdx(3,5) * shdx(2,5) * blklam
        stiff(14) = shdx(3,5) * shdy(2,5) * blklam 
        stiff(15) = shdx(3,5) * shdx(3,5) * blklam
        stiff(16) = shdx(1,5) * shdy(3,5) * blklam
        stiff(17) = shdy(1,5) * shdy(3,5) * blklam
        stiff(18) = shdx(2,5) * shdy(3,5) * blklam
        stiff(19) = shdy(2,5) * shdy(3,5) * blklam
        stiff(20) = shdx(3,5) * shdy(3,5) * blklam
        stiff(21) = shdy(3,5) * shdy(3,5) * blklam
        stiff(22) = shdx(1,5) * shdx(4,5) * blklam
        stiff(23) = shdx(4,5) * shdy(1,5) * blklam
        stiff(24) = shdx(2,5) * shdx(4,5) * blklam 
        stiff(25) = shdx(4,5) * shdy(2,5) * blklam 
        stiff(26) = shdx(3,5) * shdx(4,5) * blklam 
        stiff(27) = shdx(4,5) * shdy(3,5) * blklam 
        stiff(28) = shdx(4,5) * shdx(4,5) * blklam 
        stiff(29) = shdx(1,5) * shdy(4,5) * blklam 
        stiff(30) = shdy(1,5) * shdy(4,5) * blklam 
        stiff(31) = shdx(2,5) * shdy(4,5) * blklam 
        stiff(32) = shdy(2,5) * shdy(4,5) * blklam 
        stiff(33) = shdx(3,5) * shdy(4,5) * blklam 
        stiff(34) = shdy(3,5) * shdy(4,5) * blklam 
        stiff(35) = shdx(4,5) * shdy(4,5) * blklam 
        stiff(36) = shdy(4,5) * shdy(4,5) * blklam 
c
c..... from the viscous stiffness
c
c
      do 600 n  = 1 , 4
c$dir no_recurrence
        blkvis = evisc(n) * det(n)
ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
        stiff( 1)=stiff( 1)+( two * shdx(1,n)*shdx(1,n)
     &              +shdy (1,n) * shdy(1,n) ) *  blkvis
c
        stiff( 2)=stiff( 2)+(shdy(1,n)*shdx(1,n)) * 
     &              blkvis
c 
        stiff( 3)=stiff( 3)+( two * shdy(1,n)*shdy(1,n)
     &              +shdx(1,n)*shdx(1,n) ) * blkvis
c
        stiff( 4)=stiff( 4)+( two * shdx(1,n)*shdx(2,n)
     &              +shdy (1,n) * shdy(2,n) ) *  blkvis
c
        stiff( 5)=stiff( 5)+(shdy(2,n)*shdx(1,n)) * 
     &              blkvis
c
        stiff( 6)=stiff( 6)+( two * shdx(2,n)*shdx(2,n)
     &              +shdy(2,n) * shdy(2,n) ) * blkvis
c
        stiff( 7)=stiff( 7)+(shdy(1,n)*shdx(2,n)) * 
     &              blkvis
c
        stiff( 8)=stiff( 8)+( two * shdy(1,n)*shdy(2,n)
     &              +shdx(1,n) * shdx(2,n) ) * blkvis
c
        stiff( 9)=stiff( 9)+(shdy(2,n)*shdx(2,n)) * 
     &              blkvis
c
        stiff(10)=stiff(10)+( two * shdy(2,n)*shdy(2,n)
     &              +shdx(2,n) * shdx(2,n) ) * blkvis
c
        stiff(11)=stiff(11)+( two * shdx(1,n)*shdx(3,n)
     &              +shdy(1,n) * shdy(3,n) ) * blkvis
c
        stiff(12)=stiff(12)+(shdy(3,n)*shdx(1,n)) * 
     &              blkvis
c
        stiff(13)=stiff(13)+( two * shdx(2,n)*shdx(3,n)
     &              +shdy(2,n)*shdy(3,n) ) * blkvis
c
        stiff(14)=stiff(14)+(shdy(3,n)*shdx(2,n)) * 
     &              blkvis
c
        stiff(15)=stiff(15)+( two * shdx(3,n)*shdx(3,n)
     &              +shdy(3,n)*shdy(3,n) ) * blkvis
c
        stiff(16)=stiff(16)+(shdy(1,n)*shdx(3,n)) * 
     &              blkvis
c
        stiff(17)=stiff(17)+( two * shdy(3,n)*shdy(1,n)
     &              +shdx(3,n)*shdx(1,n) ) * blkvis
c
        stiff(18)=stiff(18)+(shdy(2,n)*shdx(3,n)) * 
     &              blkvis
c
        stiff(19)=stiff(19)+( two * shdy(3,n)*shdy(2,n)
     &              +shdx(3,n)*shdx(2,n) ) * blkvis
c
        stiff(20)=stiff(20)+(shdy(3,n)*shdx(3,n)) * 
     &              blkvis
c
        stiff(21)=stiff(21)+( two * shdy(3,n)*shdy(3,n)
     &              +shdx(3,n)*shdx(3,n) ) * blkvis
c
        stiff(22)=stiff(22)+( two * shdx(1,n)*shdx(4,n)
     &              +shdy(1,n)*shdy(4,n) ) * blkvis
c
        stiff(23)=stiff(23)+(shdy(4,n)*shdx(1,n)) * 
     &              blkvis
c
        stiff(24)=stiff(24)+( two * shdx(2,n)*shdx(4,n)
     &              +shdy(2,n)*shdy(4,n) ) * blkvis
c
        stiff(25)=stiff(25)+(shdy(4,n)*shdx(2,n)) * 
     &              blkvis
c
        stiff(26)=stiff(26)+( two * shdx(3,n)*shdx(4,n)
     &              +shdy(3,n)*shdy(4,n) ) * blkvis
c
        stiff(27)=stiff(27)+(shdy(4,n)*shdx(3,n)) * 
     &              blkvis
c
        stiff(28)=stiff(28)+( two * shdx(4,n)*shdx(4,n)
     &              +shdy(4,n)*shdy(4,n) ) * blkvis
c
        stiff(29)=stiff(29)+(shdy(1,n)*shdx(4,n)) * 
     &              blkvis
c
        stiff(30)=stiff(30)+( two * shdy(1,n)*shdy(4,n)
     &              +shdx(1,n)*shdx(4,n) ) * blkvis
c
        stiff(31)=stiff(31)+(shdy(2,n)*shdx(4,n)) * 
     &              blkvis
c
        stiff(32)=stiff(32)+( two * shdy(2,n)*shdy(4,n)
     &              +shdx(2,n)*shdx(4,n) ) * blkvis
c
        stiff(33)=stiff(33)+(shdy(3,n)*shdx(4,n)) * 
     &              blkvis
c
        stiff(34)=stiff(34)+( two * shdy(3,n)*shdy(4,n)
     &              +shdx(3,n)*shdx(4,n) ) * blkvis
c
        stiff(35)=stiff(35)+(shdy(4,n)*shdx(4,n)) * 
     &              blkvis
c
        stiff(36)=stiff(36)+( two * shdy(4,n)*shdy(4,n)
     &              +shdx(4,n)*shdx(4,n) ) * blkvis
c
c.... end block loop
c
500    continue
c
c....  end intergration loop
c
600   continue 
c
c..... adjust right hand side for a boundary condition
c
      call adjbc(lmv, stiff, ivel, vbcl, vbcr)
c 
c..... assemble this block's element stiffness
c
      call g2assm(lmv, idiag , stiff  , vlhs, ivel)
c
c.... end loop over element blocks
c
1000  continue
c
c....  return
c
      return
      end
