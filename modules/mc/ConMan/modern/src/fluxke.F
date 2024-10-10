      subroutine fluxke(x   , v   , t   , snnu)
 
c----------------------------------------------------------------------------
c
c  subroutine for time series information, this routine
c  calculates the kinetic energy and flux out the top and bottom of the box
c  at every time step using the pressure smoothing techneque.  
c  pmass is the same as tmass.
c
c----------------------------------------------------------------------------
c
c
      use fe
      use mstres
      use materials

      implicit double precision (a-h,o-z)
c
      include 'common.h'
c
      dimension x(2,*)     , v(2,*)      , t(*)     , snnu(*),
     &          tlw(4) ,    evisc(5)

      dimension rhsl(4), tmass(4)   , xl(2,4)     , vl(8)    , 
     &          prtary(8), rl(4), tl(4), tlr(4)   
      dimension vx(numnp), vz(numnp)
c
c.... initialization
c
      do i=1,numnp
        snnu(i)  = zero
        pmass(i) = zero
        vx(i) = zero
        vz(i) = zero
      enddo

      eng = zero
      tbar = zero
      ebar = zero
      total_volume = zero
      do 1000 iel = 1 , numel
c
c..... localize coordinates, temperature and  velocity
c
        do node = 1, nen
          vl(2*node-1) = v(1,ien(iel,node) ) 
          vl(2*node  ) = v(2,ien(iel,node) ) 
          xl(1,node)   = x(1,ien(iel,node) ) 
          xl(2,node)   = x(2,ien(iel,node) ) 
          tl(node)     = t(ien(iel,node)) 
          tlr(node)    = tref(ien(iel,node))
          rl(node)     = dens(ien(iel,node))
        enddo
c
c..... add adiabats to the tl
c this is needed because of the formulation in the compressible benchmark
         if (itype .lt. 3) then
          do node = 1, nen
            t_ref = T0 / diff_T
     &       * (exp(Di*(one - xl(2,node)))-one)
            tlw(node) = tl(node) + t_ref
          enddo
        else
          do node = 1, nen
            tlw(node) = tl(node) 
          enddo
        endif
c
c..... calculate kinetic energy and volume
c
        eeng= rl(1) * ( vl(1)**2 + vl(2)**2 )
     &      + rl(2) * ( vl(3)**2 + vl(4)**2 )
     &      + rl(3) * ( vl(5)**2 + vl(6)**2 )
     &      + rl(4) * ( vl(7)**2 + vl(8)**2 )
        ttbar=tlw(1) + tlw(2) + tlw(3) + tlw(4) 
c    
c.... call the global shape function
c
        volume=det(iel,1)+det(iel,2)+det(iel,3)+det(iel,4)
        flux = zero
        eng = eng + eeng * volume * pt25
        tbar= tbar+ ttbar* volume * pt25
        total_volume = total_volume + volume
c
c..... calculate flux term based on J&M(1980), eq. 52
c
        call rheol (xl, vl, tlw, iel, evisc)
        vxq = zero
        vzq = zero
        do intp=1,4
          vxq = vxq + (shl(1,intp)*vl(1) + shl(2,intp)*vl(3)
     &      + shl(3,intp)*vl(5) + shl(4,intp)*vl(7))*det(iel,intp)
          vzq = vzq + (shl(1,intp)*vl(2) + shl(2,intp)*vl(4)
     &      + shl(3,intp)*vl(6) + shl(4,intp)*vl(8))*det(iel,intp)
        enddo
        eebar=evisc(1) + evisc(2) + evisc(3) + evisc(4)
        ebar= ebar+ eebar* volume * pt25
        vxq = vxq/volume
        vzq = vzq/volume
        do intp=1, nipt-1
cc.... thermal diffusion term
          flux = flux - 
     &          (shdy(iel,1,intp)*tlw(1) + shdy(iel,2,intp)*tlw(2)
     &         + shdy(iel,3,intp)*tlw(3) + shdy(iel,4,intp)*tlw(4))
     &                 * det(iel,intp)
cc.... thermal convection term
          bdens = rl(1) * shl(1, intp) + rl(2) * shl(2,intp)
     &          + rl(3) * shl(3, intp) + rl(4) * shl(4,intp)
          flux = flux +  bdens *
     &         ( shl(1,intp)*tlw(1) + shl(2,intp)*tlw(2)
     &          +shl(3,intp)*tlw(3) + shl(4,intp)*tlw(4) )*
     &          (shl(1,intp)*vl(2) + shl(2,intp)*vl(4) +
     &           shl(3,intp)*vl(6) + shl(4,intp)*vl(8) )
     &                   *det(iel,intp)
cc.... viscous dissipation term
          slamda = - two / three
          Bvisc  = slamda + two / three
          if (itype .eq. 3) then
            Bvisc = two / three   ! compressible term is neglected for EBA
          endif
c
ccc... Tau_31 * Vz (Vx ?)  
          if (ra(mat(iel)) .gt. eps7) then
          flux = flux - (Di / ra(mat(iel))) * evisc(intp)
     &         *(shdx(iel,1,intp)*vl(2) + shdx(iel,2,intp)*vl(4)
     &         + shdx(iel,3,intp)*vl(6) + shdx(iel,4,intp)*vl(8)
     &         + shdy(iel,1,intp)*vl(1) + shdy(iel,2,intp)*vl(3)
     &         + shdy(iel,3,intp)*vl(5) + shdy(iel,4,intp)*vl(7))
     &         * (shl(1,intp)*vl(1) +  shl(2,intp)*vl(3)
     &         +  shl(3,intp)*vl(5) +  shl(4,intp)*vl(7))
     &         *  det(iel,intp)  
ccc... Tau_33 * Vz
          flux = flux - (Di / ra(mat(iel))) * evisc(intp)
     &         * (shdy(iel,1,intp)*vl(2) + shdy(iel,2,intp)*vl(4)
     &         + shdy(iel,3,intp)*vl(6) + shdy(iel,4,intp)*vl(8)
     &         + shdy(iel,1,intp)*vl(2) + shdy(iel,2,intp)*vl(4)
     &         + shdy(iel,3,intp)*vl(6) + shdy(iel,4,intp)*vl(8)
     &         + (Bvisc - two / three)
     &         *(shdx(iel,1,intp)*vl(1) + shdx(iel,2,intp)*vl(3)
     &         + shdx(iel,3,intp)*vl(5) + shdx(iel,4,intp)*vl(7)
     &         + shdy(iel,1,intp)*vl(2) + shdy(iel,2,intp)*vl(4)
     &         + shdy(iel,3,intp)*vl(6) + shdy(iel,4,intp)*vl(8)))
     &         * ( shl(1,intp)*vl(2) +  shl(2,intp)*vl(4)
     &         +  shl(3,intp)*vl(6) +  shl(4,intp)*vl(6))
     &         *  det(iel,intp)         
        endif
        enddo
c
        elnu = flux / volume
c
c smooth elnu (over elements) to snnu (over nodes)
c
        do node = 1, nen
c
c calculate element "weight"
c
          rhsl(node)=det(iel,1)*shl(node,1)+det(iel,2)*shl(node,2)
     &              +det(iel,3)*shl(node,3)+det(iel,4)*shl(node,4)
c
c assemble element stress contribution to the node
c
          snnu(ien(iel,node)) = snnu(ien(iel,node))+rhsl(node)*elnu
          vx(ien(iel,node)) = vx(ien(iel,node))+rhsl(node)*vxq
          vz(ien(iel,node)) = vz(ien(iel,node))+rhsl(node)*vzq
       enddo
c
c calculate pmass
c
        smass = rhsl(1) + rhsl(2) + rhsl(3) + rhsl(4)
        do node = 1 , nen
          tmass(node) = rhsl(node) * volume / smass
          pmass(ien(iel,node)) = pmass(ien(iel,node))+tmass(node)
        enddo
c
c end loop over elements
c
1000  continue
c
      do i=1, numnp
        snnu(i) = snnu(i) / pmass(i)
        vx(i) = vx(i) / pmass(i)
        vz(i) = vz(i) / pmass(i)
      enddo
      vrms = zero
      do i=1,numnp
        vrms = vrms + vx(i)*vx(i) + vz(i)*vz(i)
      enddo
      vrms = sqrt(vrms/float(numnp))
c
c fix boundary
c
      do i=1,nodebn
        snnu(nb(1,i)) = two * snnu(nb(1,i)) - snnu(nb(2,i))
      enddo
c
c..... save output to an array to print every 100 steps
c
      xtott = 0.0D+0
      xtotb = 0.0D+0
      tempt = 0.0D+0
      tempb = 0.0D+0
      fluxt = 0.0D+0
      fluxb = 0.0D+0
      pvelt = 0.0D+0
      pvelb = 0.0D+0
c..... flux across the bottom layer
      do node=1, nodebn/2 - 1
        xtemp = x(1,nb(1,node+1) ) - x(1,nb(1,node) )
        vtemp = v(1,nb(1,node) ) + v(1,nb(1,node+1) )
        ttemp = t( nb(1,node) ) + t(  nb(1,node+1) )
        ftemp = snnu( nb(1,node) ) + snnu(  nb(1,node+1) )
        xtotb = xtotb + xtemp
        tempb = tempb + ttemp * xtemp * pt5
        fluxb = fluxb + ftemp * xtemp * pt5
        pvelb = pvelb + vtemp * xtemp * pt5
c       write(6,*) nb(1,node), nb(1,node+1), snnu(nb(1,node)),
c    &             snnu(nb(1,node+1)), "bottom"
      enddo
c..... flux across the top layer
      do node=nodebn/2 + 1, nodebn -1
        xtemp = x(1,nb(1,node+1) ) - x(1,nb(1,node) )
        vtemp = v(1,nb(1,node) ) + v(1,nb(1,node+1) )
        ttemp = t( nb(1,node) ) + t( nb(1,node+1) )
        ftemp = snnu( nb(1,node) ) + snnu( nb(1,node+1) )
        xtott = xtott + xtemp
        tempt = tempt + ttemp * xtemp * pt5
        fluxt = fluxt + ftemp * xtemp * pt5
        pvelt = pvelt + vtemp * xtemp * pt5
c       write(6,*) nb(1,node), nb(1,node+1), snnu(nb(1,node)),
c    &             snnu(nb(1,node+1)), "top"
      enddo
      xtott = xtott+0.000001
      xtotb = xtotb+0.000001
      vrms_old = sqrt(eng)
      prtary(1) = fluxt/xtott
      prtary(2) = fluxb/xtotb
      prtary(3) = vrms
      prtary(4) = time
      prtary(5) = pvelt/xtott
      prtary(6) = pvelb/xtotb
      prtary(7) = tbar/total_volume
      prtary(8) = ebar/total_volume
      write(itsout,1001) prtary(1),prtary(2),prtary(3),
     &                   prtary(4),prtary(5),prtary(6), 
     &                   prtary(7),prtary(8),
     &                   vrms_old, lstep
c
c..... return
c
1001  format ( 9(f13.6,1x),i10)
      return
      end
