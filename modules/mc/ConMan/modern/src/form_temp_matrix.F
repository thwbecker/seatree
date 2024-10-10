      subroutine form_temp_matrix(
     &             tlhsa, tlhsb, idiag, iel)
c---------------------------------------------------------------------
c
c  This routine calculates the L.H.S. temperature matrix
c  
c----------------------------------------------------------------------
c
      use fe
      use materials
      use variables
      use sub_bench
      implicit double precision (a-h,o-z)

      include 'common.h'

c
c
      dimension  idiag(*)     , tlhsa(*)     , tlhsb(*) 
c

      common /temp2 / 
     &           el_rhs(4)    , pg(4,5)   , uq(5),   vq(5), 
     &           txq(5)       , tyq(5)       , tq(5)      ,
     &           tdq(5)       , tl(4)        , tdl(4)     ,
     &           xl(2,4)      , vl(2,4)      , uxq(5)     ,
     &           vyq(5)       , uyq(5)       , vxq(5)     ,
     &           tyqr(5)      , tlr(4)       ,  tlhsl(4,5) ,
     &           rl(4)        , evisc(5)     , blkmhu   ,
     &           conv         , tmass        , adiff    , 
     &           unorm    , eta        , xse         , blkdif   , 
     &           uxse     , ueta       , tlhsl2(4,5) , tlhsl3(4,5) ,
     &           exl(4)   , exq(5)     , tlw(4)      ,
     &           m        , k        , igadd

      do i = 1, nen
         do j = 1, nen
            tlhsl(i,j) = zero
            tlhsl2(i,j) = zero
         enddo
      enddo

c***************************************************************************
      if ( (itype .eq. 1) .or. 
     &       (itype .eq. 2) .or. 
     &       (itype .eq. 4)      ) then   ! ALA, TALA and BA

      blkmhu = dmhu(mat(iel))
c
c a hack to implement a fault in the to 50 km of the grid for
c the 2008 subduction benchmark paper
c 
      is_a_fault = 0
      if(sb_is_wedge_benchmark.eq.1)then
         call fault (v, xl, uq, vq, is_a_fault, sb_box_height,
     &        sb_plate_thickness)
      endif
      if (sb_is_wedge_benchmark .eq. 0 .or. is_a_fault.eq.0) then
c
c default case: no fault
c
        do intp = 1, nipt-1
c
c add advection with adiabatic term part of l.h.s.
c

          bdens = rl(1) * shl(1,intp) + rl(2) * shl(2,intp)
     &          + rl(3) * shl(3,intp) + rl(4) * shl(4,intp)
          do nodej = 1, nen
            conv = (shdx(iel,nodej,intp)*uq(intp) 
     &           +  shdy(iel,nodej,intp)*vq(intp)
     &           + Di*shl(nodej,intp)*vq(intp))*det(iel,intp) 

            do nodei = 1, nen
              tlhsl(nodei,nodej)  = tlhsl(nodei,nodej) + 
     &                pg(nodei,intp) * conv * bdens
            enddo
          enddo
c
c add diffusion part of l.h.s.
c
          do nodej = 1, nen
            do nodei = 1, nen
              tlhsl(nodei,nodej) = tlhsl(nodei,nodej)
     &          +(shdx(iel,nodej,intp)*shdx(iel,nodei,intp)
     &          + shdy(iel,nodej,intp)*shdy(iel,nodei,intp))
     &          *det(iel,intp)*blkdif
            enddo
          enddo
c
c end integration point loop
c
        enddo
      else
c
c this element has a fault so we have to assemble velocities carefully
c
        do intp = 1, nipt-1
c
c add advection with adiabatic term part of l.h.s.
c

           bdens = rl(1) * shl(1,intp) + rl(2) * shl(2,intp)
     &           + rl(3) * shl(3,intp) + rl(4) * shl(4,intp)
           do nodej = 1, nen
             conv = (shdx(iel,nodej,intp)*uq(intp) 
     &            +  shdy(iel,nodej,intp)*vq(intp)
     &            + Di*shl(nodej,intp)*vq(intp))*det(iel,intp) 
c
c  in the hanging wall of the fault (for node 3--above the fault-- only gauss point 1 moves)
c
             if (nodej .eq. 3) then
               conv = (shdx(iel,nodej,1)*uq(1) 
     &              +  shdy(iel,nodej,1)*vq(1)
     &              + Di*shl(nodej,1)*vq(1))*det(iel,intp) 
             endif
             do nodei = 1, nen
               tlhsl(nodei,nodej)  = tlhsl(nodei,nodej) + 
     &                 pg(nodei,intp) * conv * bdens
             enddo
           enddo
c
c add diffusion part of l.h.s.
c
          do nodej = 1, nen
            do nodei = 1, nen
              tlhsl(nodei,nodej) = tlhsl(nodei,nodej)
     &          +(shdx(iel,nodej,intp)*shdx(iel,nodei,intp)
     &          + shdy(iel,nodej,intp)*shdy(iel,nodei,intp))
     &          *det(iel,intp)*blkdif
            enddo
          enddo
c
c end integration point loop
c
         enddo
c
c end fault 
c
       endif
c
c update temperature r.h.s. residual
c
      do node = 1, nen
        el_rhs(node)  =  el_rhs(node)
     &                - tlhsl(node,1)*tl(1)
     &                - tlhsl(node,2)*tl(2)
     &                - tlhsl(node,3)*tl(3)
     &                - tlhsl(node,4)*tl(4)
      enddo

c
c update additional diffusion on r.h.s residual

      do intp = 1, nipt-1
        do nodej = 1, nen
          do nodei = 1, nen
            tlhsl2(nodei,nodej) = tlhsl2(nodei,nodej)
     &                  + shdy(iel,nodej,intp)*shdy(iel,nodei,intp)
     &                  * det(iel,intp)*blkdif
          enddo
        enddo
      enddo

      do node = 1, nen
        el_rhs(node)  =  el_rhs(node)
     &                - tlhsl2(node,1)*tlr(1)
     &                - tlhsl2(node,2)*tlr(2)
     &                - tlhsl2(node,3)*tlr(3)
     &                - tlhsl2(node,4)*tlr(4)
      enddo

c Di : dissipation number (~0.5 is typical value.)
c lamda = Kb - 2 / 3 * blkvis
c lamda: second coefficient of viscosity (hear, slamda)
c Kb    : bulk viscosity (non-dimensionalized) (hear, Bvisc)

c for Newtonian, incompressible & no dissipation (BA)
c Di = 0.0
c slamda = 0.0 (Kb = 2 / 3 * blkvis)

c for Newtonian, incompressible & dissipation (EBA)
c Di = non-zero
c slamda = 0.0 (Kb = 2 / 3 * blkvis)

c for Newtonian, compressible & dissipation (Stokesian fluid, TALA)
c Di = non-zero
c slamda = - 2 / 3 * blkvis (Kb = 0.0)
     
c for General compressible fluid
c Di = non-zero
c slamda = constant * blkvis except 0.0 and - 2 / 3 * blkvis  
c get evisc
   
      call rheol(xl, vl, tlw, iel, evisc)
            
      do intp = 1, nipt-1
c viscous dissipation term

       slamda = -two / three

       Bvisc  = slamda + two / three 

       phi = 0.0
       if (ra(mat(iel)) .gt. zero) then    
           phi = (Di/ra(mat(iel))) * evisc(intp) *
     &     ( 
     &     two * (uxq(intp) * uxq(intp) + vyq(intp) * vyq(intp)) 
     &         + (vxq(intp) + uyq(intp)) * (vxq(intp) + uyq(intp)) 
     &     + (Bvisc - two / three) 
     &     * (uxq(intp) + vyq(intp)) * (uxq(intp) + vyq(intp))
     &     )
        endif

          bdens = rl(1) * shl(1,intp) + rl(2) * shl(2,intp)
     &          + rl(3) * shl(3,intp) + rl(4) * shl(4,intp)

c for the benchmark I did not include the Di^2 term
c         tmass =  bdens*blkmhu + Di*Di*tlr(intp) + phi  
          tmass =  bdens*blkmhu + phi  

          do node=1,nen
            el_rhs(node) = el_rhs(node)
     &                   + det(iel,intp) * tmass * pg(node,intp)
        enddo

      enddo


c***************************************************************************
      else if (itype .eq. 3) then !  EBA

      blkmhu = dmhu(mat(iel))

      do intp = 1 , nipt-1
c
c add convection and adiabatic term part of l.h.s.
c
        do nodej = 1, nen
          conv = ( shdx(iel,nodej,intp)*uq(intp) 
     &         +   shdy(iel,nodej,intp)*vq(intp) ) * det(iel,intp) 
          do nodei = 1, nen
             tlhsl(nodei,nodej)  = tlhsl(nodei,nodej) + 
     &                 pg(nodei,intp) * conv 
          enddo
        enddo
c
c... adiabatic term
c
        do nodej = 1, nen
          adiab = Di*shl(nodej,intp)*vq(intp)* det(iel,intp) 
          do nodei = 1, nen
             tlhsl(nodei,nodej)  = tlhsl(nodei,nodej) + 
     &                 pg(nodei,intp) * adiab
          enddo
        enddo

c... T0/diff_T term in the adiabatic term

        do nodej = 1, nen
          adiab1 = Di*shl(nodej,intp)*vq(intp)* det(iel,intp) 
          do nodei = 1, nen
             tlhsl2(nodei,nodej)  = tlhsl2(nodei,nodej) + 
     &                 pg(nodei,intp) * adiab1
          enddo
        enddo
c      
c add diffusion part of l.h.s.
c
        do nodej = 1, nen
          do nodei = 1, nen
            tlhsl(nodei,nodej) = tlhsl(nodei,nodej)
     &        +(shdx(iel,nodej,intp)*shdx(iel,nodei,intp)
     &        + shdy(iel,nodej,intp)*shdy(iel,nodei,intp))
     &        *det(iel,intp)*blkdif
          enddo
c          el_rhs(nodej) = el_rhs(nodej)
c     &            +shl(nodej,intp)*blkmhu*det(iel,intp)
        enddo

c
c end integration point loop
c
      enddo
c
c update temperature r.h.s. residual
c
      do node = 1, nen
        el_rhs(node)  =  el_rhs(node)
     &                - tlhsl(node,1)*tl(1)
     &                - tlhsl(node,2)*tl(2)
     &                - tlhsl(node,3)*tl(3)
     &                - tlhsl(node,4)*tl(4)
      enddo

c add the T0/diff_T term to r.h.s residual

      do node = 1, nen
        el_rhs(node)  =  el_rhs(node)
     &                - tlhsl2(node,1)*tlr(1)
     &                - tlhsl2(node,2)*tlr(2)
     &                - tlhsl2(node,3)*tlr(3)
     &                - tlhsl2(node,4)*tlr(4)
      enddo


c get evisc
   
      call rheol(xl, vl, tlw, iel, evisc)
            
      do intp = 1, nipt-1
c viscous dissipation term
     
        slamda = -two / three

        Bvisc  = slamda + two / three 
               
        phi=0.0
        if (ra(mat(iel)) .gt. zero) then    
           phi = (Di/ra(mat(iel))) * evisc(intp) *
     &     ( 
     &     two * (uxq(intp) * uxq(intp) + vyq(intp) * vyq(intp)) 
     &         + (vxq(intp) + uyq(intp)) * (vxq(intp) + uyq(intp)) 
     &     )
        endif   
        
          tmass = ( - blkmhu) - phi  

          do node=1,4
            el_rhs(node) = el_rhs(node)
     &                   - det(iel,intp) * tmass * pg(node,intp)
        enddo

      enddo

      endif

c**************************************************************************

c
c add mass term
c 
c#ifdef PICARD
c#else
      if (isolve .eq. 3) then         ! Picard

      else if (isolve .eq. 2) then    ! Implicit

       do j = 1 , 4
         do i = 1 , 4
c           tlhsl(i,j)  = tlhsl(i,j)*alpha*dt
           tlhsl(i,j)  = ( tlhsl(i,j) + tlhsl2(i,j) ) *alpha*dt
         enddo
       enddo
c
       do intp = 1,4
         do nodej = 1,4
           do nodei = 1,4
c           tmass = (shl(nodej,intp) )*det(iel,intp)
           tmass = rl(nodej)*shl(nodej,intp)*det(iel,intp)
           tlhsl(nodei,nodej)=tlhsl(nodei,nodej)+pg(nodei,intp)*tmass
           enddo
         enddo
       enddo

      endif
c
c#endif
c
c assemble nonsymmetric ta and tc
c ordering for skyline method
c
        do j = 1 , 4
          do i = 1 , 4
            igadd = 0
            k = lmt(iel,j)
            if (k .ne. 0) then
              m = lmt(iel,i)
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
