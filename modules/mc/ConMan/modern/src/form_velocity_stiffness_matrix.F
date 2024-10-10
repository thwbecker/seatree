      subroutine form_velocity_stiffness_matrix(
     &                  x     , v     , t     , idiag, 
     &                  vlhsa , vlhsb , vbcr  , vrhs )


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
      use materials
      use fe
      use bc

      implicit double precision (a-h,o-z)
c
      include 'common.h'
c
      dimension  x(nsd,*)     , v(ndof,*)   , t(*)     ,
     &           idiag(*)     , vlhsa(*)    , vlhsb(*) , rl(4)    ,
     &           vbcr(*)      , vrhs(*)     , evisc(5) , tlr(4),
     &           tlw(4)  
c
      dimension stiff(8,8) , vbcl(8), tl(8), xl(2,4), vl(8)
c
c.... initialization
      do i = 1 , nsize
        vlhsa(i) = zero
        vlhsb(i) = zero
      enddo
      do i=1, neqv
         vrhs(i) = zero
         vbcr(i) = zero
      enddo
c
c.... loop over the elements
c
      do 1000 iel = 1 , numel

c.... set up the parameters

        do n = 1, nen
          xl(1,n)     = x(1, ien(iel,n) ) 
          xl(2,n)     = x(2, ien(iel,n) ) 
          vl(2*n-1)   = v(1, ien(iel,n) ) 
          vl(2*n)     = v(2, ien(iel,n) ) 
          vbcl(2*n-1) = vbc(1, ien(iel,n) )
          vbcl(2*n)   = vbc(2, ien(iel,n) )
          tl(n)       = t(  ien(iel,n) )
          tlr(n)      = tref( ien(iel,n))
          tlw(n)      = tl(n) + tlr(n) 
          rl(n)       = dens(ien(iel,n))
        enddo
        
        call rheol (xl , vl , tlw , iel , evisc)

c... form the element stiffness for this element

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    incompressible flow terms start

c....form the dynamic pressure

        blklam = alam(mat(iel))*evisc(5)*det(iel,5)*four

c....for K11, K12, K21, and K22 

        do ncol = 1, nen
          do nrow = 1, nen
            stiff(2*nrow-1, 2*ncol-1) = 
     &         shdx(iel,nrow,5)*shdx(iel,ncol,5)*blklam
            stiff(2*nrow-1, 2*ncol  ) = 
     &         shdx(iel,nrow,5)*shdy(iel,ncol,5)*blklam
             stiff(2*nrow  , 2*ncol-1) = 
     &         shdy(iel,nrow,5)*shdx(iel,ncol,5)*blklam
             stiff(2*nrow  , 2*ncol  ) = 
     &         shdy(iel,nrow,5)*shdy(iel,ncol,5)*blklam
          enddo
        enddo
           
c..... form the viscous resistance

       do n = 1, nipt-1
         blkvis = evisc(n) * det(iel,n)
         do ncol = 1, nen
           do nrow = 1, nen
            stiff(2*nrow-1, 2*ncol-1) = stiff(2*nrow-1, 2*ncol-1)
     &           + (two * shdx(iel,nrow,n) * shdx(iel,ncol,n)
     &           + shdy(iel,nrow,n) * shdy(iel,ncol,n)) * blkvis
            stiff(2*nrow-1, 2*ncol  ) = stiff(2*nrow-1, 2*ncol  )
     &           + shdy(iel,nrow,n) * shdx(iel,ncol,n)  * blkvis  
            stiff(2*nrow  , 2*ncol-1) = stiff(2*nrow  , 2*ncol-1)
     &           + shdx(iel,nrow,n) * shdy(iel,ncol,n)  * blkvis  
            stiff(2*nrow  , 2*ncol  ) = stiff(2*nrow  , 2*ncol  )
     &           + (two * shdy(iel,nrow,n) * shdy(iel,ncol,n)
     &           + shdx(iel,nrow,n) * shdx(iel,ncol,n)) * blkvis
           enddo
         enddo
       enddo
c       
c    incompressible flow terms end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c    compressible flow terms start
c
      if ( (itype .eq. 1) .or. (itype .eq. 2) ) then     ! ALA or TALA
c.... additonal terms for dynamic pressure     
        bdens = rl(1) * shl(1,5) + rl(2) * shl(2,5)
     &        + rl(3) * shl(3,5) + rl(4) * shl(4,5)
c
        rxq = rl(1) * shdx(iel,1,5) + rl(2) * shdx(iel,2,5)
     &      + rl(3) * shdx(iel,3,5) + rl(4) * shdx(iel,4,5)
c
        ryq = rl(1) * shdy(iel,1,5) + rl(2) * shdy(iel,2,5)
     &      + rl(3) * shdy(iel,3,5) + rl(4) * shdy(iel,4,5)
c            
c this is the non symmetric term
c
        do ncol = 1, nen
          do nrow = 1, nen
            stiff(2*nrow-1, 2*ncol-1) = stiff(2*nrow-1, 2*ncol-1)
     &       + shdx(iel,nrow,5)*rxq*shl(ncol,5)*blklam/bdens
            stiff(2*nrow-1, 2*ncol  ) = stiff(2*nrow-1, 2*ncol  )
     &       + shdx(iel,nrow,5)*ryq*shl(ncol,5)*blklam/bdens
            stiff(2*nrow  , 2*ncol-1) = stiff(2*nrow  , 2*ncol-1)
     &       + shdy(iel,nrow,5)*rxq*shl(ncol,5)*blklam/bdens
            stiff(2*nrow  , 2*ncol  ) = stiff(2*nrow  , 2*ncol  )
     &       + shdy(iel,nrow,5)*ryq*shl(ncol,5)*blklam/bdens
         enddo
        enddo 
c
c.... additional terms for viscous resistance  
c
       do inpt=1,nipt-1
         blkvis = evisc(nipt) * det(iel,nipt)
         Bvisc = zero * blkvis
         slamda = Bvisc  - two / three * blkvis
         do ncol = 1, nen
          do nrow = 1, nen
            stiff(2*nrow-1, 2*ncol-1) = stiff(2*nrow-1, 2*ncol-1)
     &        + shdx(iel,nrow,inpt) * shdx(iel,ncol,inpt) * slamda
            stiff(2*nrow-1, 2*ncol  ) = stiff(2*nrow-1, 2*ncol  )
     &        + shdx(iel,nrow,inpt) * shdy(iel,ncol,inpt) * slamda
            stiff(2*nrow  , 2*ncol-1) = stiff(2*nrow  , 2*ncol-1)
     &        + shdy(iel,nrow,inpt) * shdx(iel,ncol,inpt) * slamda
            stiff(2*nrow  , 2*ncol  ) = stiff(2*nrow  , 2*ncol  )
     &        + shdy(iel,nrow,inpt) * shdy(iel,ncol,inpt) * slamda
          enddo
         enddo
        enddo        
c
c    compressible flow terms end
c
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c     
c..... adjust right hand side for a boundary condition
c  
      call adjust_velocity_rhs_for_bc
     &  (lmv, stiff, iel, vbcl, vbcr)
c 
c..... assemble this blocks element stiffness
c
      call velo_stiff_matrix_assemble
     &  (lmv, idiag , stiff , vlhsa, vlhsb, iel)
c
c.... end loop over element blocks
c
1000  continue
c
c....  return
c
      return
      end
