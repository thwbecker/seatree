       subroutine print_compbm_data(x, v, t)
c
c.... program to print data for comparison with ET, PvK and J&M
c.... temperature table corresponding to nodal values
c.... R0, Nu #, vrms, vsurfmax, <vsurf>, <T>, <phi>, <W>, 
c.... dV(%), E=<phi>/F', <E> and time                  
c
      use materials
      use fe

      implicit double precision (a-h,o-z)
c
      include 'common.h'

      dimension x(ndof,*) , v(ndof,*) , t(*)

      dimension rl(4)    ,  rhsl(4)   , tmass(4) ,
     &          uq(5)     , vq(5)     , tq(5)    , tl(4)    ,
     &          xl(2,4)   , vl(2,4)   , uxq(5)   , vxq(5)   ,
     &          vyq(5)    , uyq(5)    , evisc(5) , phi(5)   ,
     &          bdens(5)  , 
     &          tlr(4)   , tlw(4)  


      dimension pmass(numnp), vx(numnp), vz(numnp)
 
c.... initialization

      do i = 1, numnp
        pmass(i) = zero; vx(i) = zero ; vz(i) = zero
      enddo
      do i = 1, 4
        rl(i) = zero
        tmass(i) = zero
        rhsl(i) = zero
        phi(i) = zero
      enddo
   
c     nltop = nodal number of left-top edge
c     ndif = difference of nodal number between adjacent column    

      nrow   = nelz + 1
      ncol   = nelx + 1
      pflux  = one
      elen_x = x(1,numnp) / dble(nelx)
      elen_z = x(2,numnp) / dble(nelz)
      R0     = ra(mat(1))
      nltop  = nrow
      ndif   = nrow    

c      e_vol  = elen_x * elen_z  -> calculated by Jacobian

c.... 1. display input data
   
c.... find vsurfmax among x-velocity values for the top surface
c.... vsurfmax is absolute value
         
      vsurfm = dabs(v(1,nltop))
      do i = nltop,numnp,ndif
         vtemp  = dabs(v(1,i))
         if (vsurfm .lt. vtemp) then
           vsurfm = vtemp
         endif
      enddo        
      
c.... find <vsurf> using trapezoid rule
      
      vsurf = zero
      do i  = nltop,numnp,ndif
        vsurf = vsurf + v(1,i)
      enddo
      vsurf = vsurf - pt5 * (v(1,nltop) + v(1,numnp))
      vsurf = vsurf * elen_x

c.... find <vrms> for root mean square of velocity       

      do iel = 1, numel

       do n=1,nen
        xl(1,n) = x( 1,ien(iel,n) ) 
        xl(2,n) = x( 2,ien(iel,n) ) 
        vl(1,n) = v( 1,ien(iel,n) ) 
        vl(2,n) = v( 2,ien(iel,n) )
       enddo

       vlx = zero
       vlz = zero

       do intp = 1, nipt-1
          vlx  = vlx + ( vl(1,1)*shl(1,intp) + vl(1,2)*shl(2,intp)
     &                 + vl(1,3)*shl(3,intp) + vl(1,4)*shl(4,intp) )
     &                 * det(iel,intp)
          vlz  = vlz + ( vl(2,1)*shl(1,intp) + vl(2,2)*shl(2,intp)
     &                 + vl(2,3)*shl(3,intp) + vl(2,4)*shl(4,intp) )
     &                 * det(iel,intp)
        enddo

        volume = det(iel,1)   + det(iel,2)   + det(iel,3)   + det(iel,4)
        vlx = vlx / volume
        vlz = vlz / volume

c calculate element "weight"

        do n = 1, nen
          rhsl(n) = det(iel,1) * shl(n,1) + det(iel,2) * shl(n,2)
     &            + det(iel,3) * shl(n,3) + det(iel,4) * shl(n,4)
          vx(ien(iel,n)) = vx(ien(iel,n)) + rhsl(n) * vlx
          vz(ien(iel,n)) = vz(ien(iel,n)) + rhsl(n) * vlz
        enddo

        do n = 1, nen
          tmass(n) = det(iel,1) * shl(n,1) + det(iel,2) * shl(n,2)
     &             + det(iel,3) * shl(n,3) + det(iel,4) * shl(n,4)
        enddo

        smass = tmass(1) + tmass(2) + tmass(3) + tmass(4)

        do n = 1, nen
          tmass(n) = tmass(n) * volume / smass
        enddo

        do n = 1, nen
          pmass(ien(iel,n)) = pmass(ien(iel,n)) + tmass(n)
        enddo
 
       enddo

        do i=1, numnp
          vx(i) = vx(i) / pmass(i)
          vz(i) = vz(i) / pmass(i)
        enddo

      vrms = zero
      do i = 1, numnp
        vrms = vrms + ( vx(i) * vx(i) + vz(i) * vz(i) )
      enddo
      vrms = dsqrt(vrms / numnp) 

c.... find <T> using FEM volume integral
      
      T_mean = zero
      p_mean = zero
      w_mean = zero
      H_mean = zero
      E_mean = zero

c.... sum of phi, heat generation, work calculated from each element
c.... order is from the left bottom to the left top by moving from left to right
c.... the same of nodal ordering

      do iel = 1, numel        

        do n = 1 , nen
          xl(1,n) = x(1, ien(iel,n))
          xl(2,n) = x(2, ien(iel,n))
          vl(1,n) = v(1, ien(iel,n))
          vl(2,n) = v(2, ien(iel,n))
          rl(n)   = dens(ien(iel,n))
          tl(n)   = t(   ien(iel,n))
          tlr(n)  = tref(ien(iel,n))
          tlw(n)  = tl(n) + tlr(n) 
        enddo

c....   calculation of element volume
        e_vol = det(iel,1) + det(iel,2) + det(iel,3) + det(iel,4)

c....   calculation of density 
        do intp = 1, nipt
          bdens(intp) = rl(1) * shl(1,intp) + rl(2) * shl(2,intp)
     &                + rl(3) * shl(3,intp) + rl(4) * shl(4,intp)
        enddo

c....   generate related values

        do intp = 1, nipt
          uq(intp)  = vl(1,1) * shl(1,intp) + vl(1,2) * shl(2,intp) 
     &              + vl(1,3) * shl(3,intp) + vl(1,4) * shl(4,intp)

          vq(intp)  = vl(2,1) * shl(1,intp) + vl(2,2) * shl(2,intp) 
     &              + vl(2,3) * shl(3,intp) + vl(2,4) * shl(4,intp)

          tq(intp)  = tlw(1) * shl(1,intp) + tlw(2) * shl(2,intp)
     &              + tlw(3) * shl(3,intp) + tlw(4) * shl(4,intp)

          uxq(intp)=vl(1,1)*shdx(iel,1,intp)+vl(1,2)*shdx(iel,2,intp)
     &             +vl(1,3)*shdx(iel,3,intp)+vl(1,4)*shdx(iel,4,intp)

          vyq(intp)=vl(2,1)*shdy(iel,1,intp)+vl(2,2)*shdy(iel,2,intp)
     &             +vl(2,3)*shdy(iel,3,intp)+vl(2,4)*shdy(iel,4,intp)

          uyq(intp)=vl(1,1)*shdy(iel,1,intp)+vl(1,2)*shdy(iel,2,intp)
     &             +vl(1,3)*shdy(iel,3,intp)+vl(1,4)*shdy(iel,4,intp)

          vxq(intp)=vl(2,1)*shdx(iel,1,intp)+vl(2,2)*shdx(iel,2,intp)
     &             +vl(2,3)*shdx(iel,3,intp)+vl(2,4)*shdx(iel,4,intp)

        enddo

c....   calculation of viscous dissipation      
c....   thermal expansivity variation, exq is divided by the same
c....   from Ra number.

        call rheol(xl, vl, tlw, iel, evisc)

        slamda = -two / three
        Bvisc  = slamda + two / three 
        if ( (itype .eq. 1) .or. (itype .eq. 2) ) then   ! ALA and TALA

        do intp = 1, nipt
          phi(intp) = (Di / ra(mat(iel))) * evisc(intp) * 
     &    (two * (uxq(intp) * uxq(intp) + vyq(intp) * vyq(intp)) 
     &         + (vxq(intp) + uyq(intp))*(vxq(intp) + uyq(intp)) 
     &    + (Bvisc - two / three) 
     &    * (uxq(intp) + vyq(intp)) * (uxq(intp) + vyq(intp)))        
        enddo

        else if (itype .eq. 3) then                      ! EBA

        if (ra(mat(iel)) .gt. eps7) then
        do intp = 1, nipt   
           phi(intp) = (Di / ra(mat(iel))) * evisc(intp) * 
     &    (two * (uxq(intp) * uxq(intp) + vyq(intp) * vyq(intp)) 
     &         + (vxq(intp) + uyq(intp)) * (vxq(intp) + uyq(intp)))          
        enddo
        endif
        endif

c....   calculation of internal heat generation

        blkmhu = dmhu(mat(iel)) 

c....   calculation of each value for each element

        T_elem = zero
        phi_el = zero
        w_el   = zero
        H_elem = zero
        E_elem = zero

        do n = 1, nen
          T_elem = T_elem + (tq(n) - T0/diff_T) * shl(n,5)
          phi_el = phi_el + phi(n)  * shl(n,5)
          w_el   = w_el   +  Di * bdens(n)
     &           * tq(n) * vq(n) * shl(n,5)

          H_elem = H_elem + bdens(n)* blkmhu  * shl(n,5)
          E_elem = E_elem + bdens(n)* tq(n)* shl(n,5)
        enddo

        T_mean = T_mean + T_elem * e_vol
        p_mean = p_mean + phi_el * e_vol
        w_mean = w_mean + w_el   * e_vol
        H_mean = H_mean + H_elem * e_vol 
        E_mean = E_mean + E_elem * e_vol

c.... print phi, W and viscosity corresponding to each element
c.... save in the viscosity.bm1 file (file is opened and closed in the conman.c)        
        if (mod(time, tsave) .lt. dt) then
          x_cent = (xl(1,1)+xl(1,2)+xl(1,3)+xl(1,4))/four
          y_cent = (xl(2,1)+xl(2,2)+xl(2,3)+xl(2,4))/four

          work_e = w_el   * e_vol
          phi_e  = phi_el * e_vol

          visc_e = ( evisc(1) + evisc(2) + evisc(3) + evisc(4) ) / four
        endif
      enddo

c.... calculation of dV(%) term

      if ( abs(Di) .lt. eps7 ) then
        dv = zero
      else 
        dv = (p_mean-w_mean)/w_mean*100.000000D+0
      endif

c.... calculation of Nu #, for ET and ConMan, mean T of whole domain is 0.5

       if (T_mean .gt. eps7 ) pNu = pt5 / T_mean

c.... efficiency, E = <phi> / F'. For PvK's experiments, F' = diff_T
c      effic = p_mean / diff_T
       effic = p_mean

c....  display result data
      write(icomp, 1110) R0, Di, vrms, vsurfm, vsurf, 
     &                   T_mean, p_mean, w_mean
 1110 format(8(1pe14.7,1x))

      return
c
      end
