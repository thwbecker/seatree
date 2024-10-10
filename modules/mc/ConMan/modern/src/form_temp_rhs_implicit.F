         subroutine form_temp_rhs_implicit( t     , tdot  ,
     &                        x     , v     , idiagt,
     &                        trhs  , tlhsa , tlhsb )
c----------------------------------------------------------------------
c
c  S.D.K. 03/25/92 IMPLICIT and PICARD version
c  C.L.   02/28/07 TALA added
c 
c----------------------------------------------------------------------
c
      use fe
      use bc
      use materials
      use sub_bench
      implicit double precision (a-h,o-z)
c
      include 'common.h'
c
      dimension  t(*)         , tdot(*)      ,
     &           x(ndof,*)    , v(ndof,*)    , idiagt(*)  ,
     &           trhs(*)      , tlhsa(*)     , tlhsb(*)   ,
     &           eval(6)      
        
      common /temp2 / 
     &           el_rhs(4)    , pg(4,5)   , uq(5),   vq(5),
     &           txq(5)       , tyq(5)       , tq(5)      ,
     &           tdq(5)       , tl(4)        , tdl(4)     ,
     &           xl(2,4)      , vl(2,4)      , uxq(5)     ,
     &           vyq(5)       , uyq(5)       , vxq(5)     ,
     &           tyqr(5)      , tlr(4)       , tlhsl(4,5) ,
     &           rl(4)        , evisc(5)     , blkmhu    ,
     &           conv         , tmass        , adiff     ,
     &           unorm    , eta        , xse         , blkdif   , 
     &           uxse     , ueta       , tlhsl2(4,5) , tlhsl3(4,5) ,
     &           exl(4)   , exq(5)     , tlw(4)      ,
     &           m        , k        , igadd

c.... local variables used here

c     blkmhu   ,  conv      , tmass    , adiff    , 
c     unorm    , eta        , xse      , blkdif   , 
c     uxse     , ueta       , m        , k        ,
c     igadd  

c.... initialization
      do i = 1, 4
        do j = 1, 5
          tlhsl(i,j)= zero
        enddo
      enddo

      do i = 1 , numnp
        trhs(i) = zero
      enddo

      do i = 1 , nsizet
        tlhsa(i) = zero
        tlhsb(i) = zero
      enddo
c
c.... loop over the elements
c
      do iel = 1 , numel
c
c..... localize the temperature, tdot, velocity, density
c
      do n = 1,4
        xl(1,n) = x(1, ien(iel,n) ) 
        xl(2,n) = x(2, ien(iel,n) ) 
        vl(1,n) = v(1, ien(iel,n) ) 
        vl(2,n) = v(2, ien(iel,n) )
        rl(n)   = dens(ien(iel,n) ) 
      enddo
c
      if (isolve .eq. 3) then       ! Picard
        do n = 1, 4
          tl(n) = tbc( ien(iel,n) )
          tlr(n) = tref(ien(iel,n))
          tlw(n) = t(ien(iel,n)) + tlr(n)
        enddo
      else if (isolve .eq. 2) then  ! Implicit
        do n = 1, 4
          tl(n) = t( ien(iel,n) )
          tlr(n) = tref(ien(iel,n))
          tlw(n) = t(ien(iel,n)) + tlr(n)
        enddo
      endif
c
c....form the the values of the temperature and velocity at intergration points
c
        do intp = 1,5
          uq(intp)  = vl(1,1) * shl(1,intp) + vl(1,2) * shl(2,intp) 
     &              + vl(1,3) * shl(3,intp) + vl(1,4) * shl(4,intp)

          vq(intp)  = vl(2,1) * shl(1,intp) + vl(2,2) * shl(2,intp) 
     &              + vl(2,3) * shl(3,intp) + vl(2,4) * shl(4,intp) 

          tq(intp)  =tl(1) * shl(1,intp) + tl(2) * shl(2,intp)
     &              +tl(3) * shl(3,intp) + tl(4) * shl(4,intp)

          txq(intp) =tl(1)*shdx(iel,1,intp)+tl(2)*shdx(iel,2,intp)
     &              +tl(3)*shdx(iel,3,intp)+tl(4)*shdx(iel,4,intp)

          tyq(intp) =tl(1)*shdy(iel,1,intp)+tl(2)*shdy(iel,2,intp)
     &              +tl(3)*shdy(iel,3,intp)+tl(4)*shdy(iel,4,intp)

          uxq(intp) =vl(1,1)*shdx(iel,1,intp)+vl(1,2)*shdx(iel,2,intp)
     &              +vl(1,3)*shdx(iel,3,intp)+vl(1,4)*shdx(iel,4,intp)

          vyq(intp) =vl(2,1)*shdy(iel,1,intp)+vl(2,2)*shdy(iel,2,intp)
     &              +vl(2,3)*shdy(iel,3,intp)+vl(2,4)*shdy(iel,4,intp)

          uyq(intp) =vl(1,1)*shdy(iel,1,intp)+vl(1,2)*shdy(iel,2,intp)
     &              +vl(1,3)*shdy(iel,3,intp)+vl(1,4)*shdy(iel,4,intp)

          vxq(intp) =vl(2,1)*shdx(iel,1,intp)+vl(2,2)*shdx(iel,2,intp)
     &              +vl(2,3)*shdx(iel,3,intp)+vl(2,4)*shdx(iel,4,intp)
 
          tyqr(intp)=tlr(1) *shdy(iel,1,intp)+tlr(2) *shdy(iel,2,intp)
     &              +tlr(3) *shdy(iel,3,intp)+tlr(4) *shdy(iel,4,intp) 
        enddo
        is_a_fault = 0
        if(sb_is_wedge_benchmark.eq.1)then
           call fault (v, xl, uq, vq, is_a_fault, sb_box_height,
     &          sb_plate_thickness)
        endif
c
c...for PG shape functions (see Brooks Thesis)
c
c.... calculate min element dimension per node
c
        exse1 = pt5*( xl(1,2) + xl(1,3) - xl(1,4) - xl(1,1) )
        exse2 = pt5*( xl(2,2) + xl(2,3) - xl(2,4) - xl(2,1) )
        eeta1 = pt5*( xl(1,3) + xl(1,4) - xl(1,1) - xl(1,2) )
        eeta2 = pt5*( xl(2,3) + xl(2,4) - xl(2,1) - xl(2,2) )
        hxse = dsqrt(exse1*exse1 + exse2*exse2)
        heta = dsqrt(eeta1*eeta1 + eeta2*eeta2)
        eval(1) = exse1/hxse
        eval(2) = exse2/hxse
        eval(3) = eeta1/heta
        eval(4) = eeta2/heta
        eval(5) = hxse
        eval(6) = heta
        uxse = dabs(uq(5)*eval(1) + vq(5)*eval(2)) * eval(5)  
        ueta = dabs(uq(5)*eval(3) + vq(5)*eval(4)) * eval(6) 

        blkmhu = dmhu(mat(iel)) 
        blkdif = diff(mat(iel))

        if( uxse .gt. two*blkdif)  then
           xse = one - two*blkdif/uxse
        else
           xse = zero
        end if
c
        if( ueta .gt. two*blkdif)  then
           eta = one - two*blkdif/ueta
        else
           eta = zero
        end if
c
        unorm = uq(5) * uq(5) + vq(5) * vq(5)
        if( unorm .gt. 0.000001 ) then
          adiff = (uxse * xse + ueta * eta)/(two*unorm)
        else
          adiff = zero
        end if
c
c loop over integration points 
c
        do node = 1,4
          do intp = 1,5
c
c..... form pg shape function's for intergration point intp
c
           pg(node,intp) = shl(node,intp) + 
     &       adiff * (uq(intp)*shdx(iel,node,intp)
     &             +  vq(intp)*shdy(iel,node,intp))
          enddo
          el_rhs(node) = zero
        enddo
c
c....  form Residue term at intergation point INTP (done in f_tlhs)
c

         call form_temp_matrix(
     &          tlhsa, tlhsb, idiagt, iel, ifault)
c
c..... assemble this element's residual
c
        do n = 1,4
          if(lmt(iel,n) .ne. 0) 
     &    trhs(lmt(iel,n)) = trhs(lmt(iel,n)) + el_rhs(n)
        enddo
c      
c.... end loop over elements
c
      enddo

c.... adjust for heat flux boundary conditions.
c.... assumes user inputs heat-flux integrated across the boundary length 
c.... heat flux is applied on the element length unit
c.... For isuf(n,2), when n is 1:bottom, 2:right, 3:top and 4:left
c.... shl(n,5) = 0.25

      do n = 1, numsuf
        iel = isuf(n,1)
        if (isuf(n,2) .eq. 1) then
          trhs(lmt(iel,1))=trhs(lmt(iel,1))-suf(n,3)*two*shl(1,5)
          trhs(lmt(iel,2))=trhs(lmt(iel,2))-suf(n,3)*two*shl(2,5)
        end if
        if (isuf(n,2) .eq. 2) then
          trhs(lmt(iel,2))=trhs(lmt(iel,2))-suf(n,3)*two*shl(2,5)
          trhs(lmt(iel,3))=trhs(lmt(iel,3))-suf(n,3)*two*shl(3,5)
        end if
        if (isuf(n,2) .eq. 3) then
          trhs(lmt(iel,3))=trhs(lmt(iel,3))-suf(n,3)*two*shl(3,5)
          trhs(lmt(iel,4))=trhs(lmt(iel,4))-suf(n,3)*two*shl(4,5)
        end if
        if (isuf(n,2) .eq. 4) then
          trhs(lmt(iel,4))=trhs(lmt(iel,4))-suf(n,3)*two*shl(4,5)
          trhs(lmt(iel,1))=trhs(lmt(iel,1))-suf(n,3)*two*shl(1,5)
        end if
      enddo
c
      call unfact(tlhsa, tlhsb, idiagt, neqt)
      call unback(tlhsa, tlhsb, trhs, idiagt, neqt)
c
c.... adjust group assembled residual (trhs) for boundray conditions
c.... (i.e replace with boundray value)
c

      if (isolve .eq. 2) then   ! Implicit
        do n = 1 , numnp
          if(idt(n) .ne. 0 ) then
            tdot(n) = trhs(idt(n))
            t(n)    = t(n) + alpha * dt * tdot(n) 
          else
            t(n) = tbc(n)
            tdot(n) = zero
          end if
        enddo
      else if (isolve .eq. 3) then   ! Picard
        do n = 1 , numnp
          if(idt(n) .ne. 0 ) then 
c           t(n) = (t(n)+ trhs(idt(n)))/2.0e0      ! original
            t(n) = 0.7d0*t(n)+ 0.3d0*trhs(idt(n))  ! slower but, safer (avoid overshoots)
            tdot(n) = zero
          else
            t(n) = tbc(n)
            tdot(n) = zero
          end if
        enddo
      endif
c
c... wrap b.c.
c
      if (nwrap .gt. 0) then
        do n = 1, nwrap
          iwrap = numnp - nwrap + n
          t(iwrap) = t(n)
          tdot(iwrap) = tdot(n)
        enddo
      endif
c
c....  return
c
      return
      end
