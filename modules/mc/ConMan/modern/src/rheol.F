      subroutine rheol( xl , vl , tl , iel , evisc)  
c
c----------------------------------------------------------------------
c  This routine calculates the viscosities per each element unit 
c  (Changyeol Lee Sept. 06)
c  
c     input tl is non-dimensonalized with a temperature unit of Kelvin (K)
c
c     tl = t-non-dim + T0/diff_T (so we have to take that factor out for non-dim...)
c      
c----------------------------------------------------------------------
c
      use materials
      use fe

      implicit double precision (a-h,o-z)
c
      include 'common.h'
c
      dimension xl(2,4), vl(2,4), tl(4), evisc(5), tll(4)
      dimension evisc1(5), evisc2(5) 
      
      vismin = 1.0d-5           ! viscosity minimum cut
      vismax = 1.0d5            ! viscosity maximum cut

      Tfac = T0 / diff_T
      
      do i = 1, 4
       tll(i) = zero
      enddo
      
c********************************************************************
c..... constant viscosity - Newtonian fluid
c
      if ( ntimvs .eq. 0 ) then
        do intp = 1, nipt
          evisc(intp) = visc(mat(iel))
        enddo
c
c********************************************************************
c      temperature dependent rheology (Blackenbach et al. 1989)
c      case2a 
c      
      elseif ( ntimvs .eq. 1 ) then 
c                                 
c..... add adiabatic temperature profile to tl 
c
        do i = 1, 4
          tll(i) = tl(i) - Tfac ! for Blankenbach scale   
        enddo
c
        do intp = 1, nipt 
          tq = tll(1) * shl(1,intp) + tll(2) * shl(2,intp)
     &       + tll(3) * shl(3,intp) + tll(4) * shl(4,intp)     
          zq = xl(2,1) * shl(1,intp) + xl(2,2) * shl(2,intp)
     &       + xl(2,3) * shl(3,intp) + xl(2,4) * shl(4,intp)
          xq = xl(1,1) * shl(1,intp) + xl(1,2) * shl(2,intp)
     &       + xl(1,3) * shl(3,intp) + xl(1,4) * shl(4,intp)
 
c.... calculation of exponent term for creep laws
          bconst = 6.907755279d0     ! for case2a 
          cconst = 0.000000000d0
c
c.... viscosity
c
          evisc(intp) = visc(mat(iel)) 
     &                * exp( -bconst*tq+cconst*(1.0d0-zq) )
        enddo
c
c********************************************************************
c      temperature dependent rheology (Blackenbach et al. 1989)
c      case2b 
c
      elseif ( ntimvs .eq. 2 ) then
c
c..... add adiabatic temperature profile to tl 
c
        do i = 1, 4
          tll(i) = tl(i) - Tfac ! for Blankenbach scale   
        enddo
c
        do intp = 1, nipt
          tq = tll(1) * shl(1,intp) + tll(2) * shl(2,intp)
     &       + tll(3) * shl(3,intp) + tll(4) * shl(4,intp)     
          zq = xl(2,1) * shl(1,intp) + xl(2,2) * shl(2,intp)
     &       + xl(2,3) * shl(3,intp) + xl(2,4) * shl(4,intp)
          xq = xl(1,1) * shl(1,intp) + xl(1,2) * shl(2,intp)
     &       + xl(1,3) * shl(3,intp) + xl(1,4) * shl(4,intp)
 
c.... calculation of exponent term for creep laws
          bconst = 9.704060528d0     ! for case2b 
          cconst = 4.158883083d0
c
c.... viscosity
c
          evisc(intp) = visc(mat(iel)) 
     &                * exp( -bconst*tq+cconst*(1.0d0-zq) )
        enddo
c
c********************************************************************
c      temperature dependent rheology 
c      Arrhenius law (diffusion creep)
c      
      elseif ( ntimvs .eq. 3 ) then
c                                 
c..... add adiabatic temperature profile to tl 
c
         do i = 1, 4
           tll(i) = tl(i) - Tfac ! for Blankenbach scale   
         enddo
c 
        do intp = 1, nipt   
          tq = tll(1) * shl(1,intp) + tll(2) * shl(2,intp)
     &       + tll(3) * shl(3,intp) + tll(4) * shl(4,intp)     
          zq = xl(2,1) * shl(1,intp) + xl(2,2) * shl(2,intp)
     &       + xl(2,3) * shl(3,intp) + xl(2,4) * shl(4,intp)
          xq = xl(1,1) * shl(1,intp) + xl(1,2) * shl(2,intp)
     &       + xl(1,3) * shl(3,intp) + xl(1,4) * shl(4,intp)
c
c.... viscosity
c
          power = visexp(tq, tcon(1,mat(iel)),
     &                   zq, tcon(2,mat(iel)) )
c
          power0 = visexp(one, tcon(1,mat(iel)),
     &                   zero, tcon(2,mat(iel)) )
          evisc(intp) = visc(mat(iel)) * (exp(power - power0))
        enddo
      elseif (ntimvs .eq. 4) then

c      
c********************************************************************      
c..... combination of diffusion and dislocation creeps with
c      temperature and pressure dependent rheology

c.... strain rate conversion paramters   
       pkappa = 1.0d-6   ! thermal diffusivity 
       bottom = 2.8900d6  ! depth of domain
c.... constants for creeps
       ! Karato and Wu, 1993 (wet Olivine)
c       p_grain = 2.5d0    ! grain size exponent
c       grain_size  = 1.0d-3 ! grain size
c      a_dis  =  3.9 * 1.0d-15  ! constant of dislocation creep
c       a_dif   =  3.7 * 1.0d-19  ! constant of diffusion creep
c       xn = 3.0d0   ! stress exponent for dislocation creep

       ! Karato and Wu, 1993 (dry Olivine)
       p_grain = 2.5d0    ! grain size exponent
       grain_size  = 1.0d-3 ! grain size
       a_dis  =  2.4 * 1.0d-16  ! constant of dislocation creep
       a_dif   =  6.1 * 1.0d-19  ! constant of diffusion creep
       xn = 3.5d0   ! stress exponent for dislocation creep

       visrat = one / ( a_dis ** (one/xn) ) 
     &             / ( (grain_size ** p_grain) / a_dif )  ! viscosity ratio

       rate_conv = pkappa / (bottom * bottom)  ! strain rate conversion factor   6.25e-18
c
c..... add adiabatic temperature profile to tl 
c
        do i = 1, nipt
          tll(i) = tl(i) - Tfac 
        enddo

        do intp = 1, nipt   
          tq = tll(1) * shl(1,intp) + tll(2) * shl(2,intp)
     &       + tll(3) * shl(3,intp) + tll(4) * shl(4,intp)
          zq = xl(2,1) * shl(1,intp) + xl(2,2) * shl(2,intp)
     &       + xl(2,3) * shl(3,intp) + xl(2,4) * shl(4,intp)
          xq = xl(1,1) * shl(1,intp) + xl(1,2) * shl(2,intp)
     &       + xl(1,3) * shl(3,intp) + xl(1,4) * shl(4,intp)

c.... calculation of exponent term for diffusion creep 

c         call ve1(tq, tcon(1,mat(iel)),
c    &            zq, tcon(2,mat(iel)), power1)

          estar = 300.0
          vstar = 6.0
          power = visexp(tq, estar,
     &                   zq, vstar )
          power0 = visexp(one, estar,
     &                   zero, vstar )
          evisc1(intp) = visc(mat(iel)) * (exp(power - power0))
          if (evisc1(intp) .gt. 1e6) evisc1(intp) = 1.0e6

c.... calculation of exponent term for dislocation creep 

c.... find velocity gradient at element centers 

c  normal gradient toward x (Vx,x)
          edotxx = 
     &            (shdx(iel,1,5)*vl(1,1) + shdx(iel,2,5)*vl(1,2)
     &           + shdx(iel,3,5)*vl(1,3) + shdx(iel,4,5)*vl(1,4) )

c  normal gradient toward y (Vy,y)
          edotyy = 
     &            (shdy(iel,1,5)*vl(2,1) + shdy(iel,2,5)*vl(2,2)
     &           + shdy(iel,3,5)*vl(2,3) + shdy(iel,4,5)*vl(2,4) )

c  shear gradient toward x or y (Vx,y = Vy,x)
          edotxy = pt5 * 
     &            (shdx(iel,1,5)*vl(2,1) + shdx(iel,2,5)*vl(2,2)
     &           + shdx(iel,3,5)*vl(2,3) + shdx(iel,4,5)*vl(2,4)
     &           + shdy(iel,1,5)*vl(1,1) + shdy(iel,2,5)*vl(1,2)
     &           + shdy(iel,3,5)*vl(1,3) + shdy(iel,4,5)*vl(1,4))

c.... second invariant of the strain rate

          edot2 = edotxx * edotxx + edotyy * edotyy
     &            +  two * edotxy * edotxy

          secinv = dsqrt(pt5 * edot2)

          secinv = rate_conv * secinv * 1.0e-1 +1.0e-17  ! conversion edot to real value (/s) 

          evisc2(intp) = evisc1(intp) *
     &                * ( secinv ** ( (one - xn) / xn ) )

c.... total viscosity (1/diff_creep + 1/disloc_creep)^(-1)

          evisc(intp) = one/(one/evisc1(intp) + one/evisc2(intp)) 

c.... set the viscous cut 
c vismax = 1e5, vismin = 1e-5

         if (evisc(intp) .gt. vismax) then
           evisc(intp) = vismax
         endif

         if (evisc(intp) .lt. vismin) then
           evisc(intp) = vismin
         endif

       enddo
      elseif (ntimvs .eq. 5) then
c********************************************************************      
c  Stein and Hansen Rheology
c
c..... add adiabatic temperature profile to tl
c
        do i = 1, nipt
         tll(i) = tl(i) - Tfac
        enddo

        do intp = 1, nipt
          tq = tll(1) * shl(1,intp) + tll(2) * shl(2,intp)
     &       + tll(3) * shl(3,intp) + tll(4) * shl(4,intp)
          zq = xl(2,1) * shl(1,intp) + xl(2,2) * shl(2,intp)
     &       + xl(2,3) * shl(3,intp) + xl(2,4) * shl(4,intp)
          xq = xl(1,1) * shl(1,intp) + xl(1,2) * shl(2,intp)
     &       + xl(1,3) * shl(3,intp) + xl(1,4) * shl(4,intp)

c this is natural log of 10^5
          const1 = 11.512925
          evisc1(intp) = exp(-const1*tq+tcon(1,mat(iel))*(one-zq))

c.... find velocity gradient at element centers 

c  normal gradient toward x (Vx,x)
          edotxx = 
     &            (shdx(iel,1,5)*vl(1,1) + shdx(iel,2,5)*vl(1,2)
     &           + shdx(iel,3,5)*vl(1,3) + shdx(iel,4,5)*vl(1,4) )

c  normal gradient toward y (Vy,y)
          edotyy = 
     &            (shdy(iel,1,5)*vl(2,1) + shdy(iel,2,5)*vl(2,2)
     &           + shdy(iel,3,5)*vl(2,3) + shdy(iel,4,5)*vl(2,4) )

c  shear gradient toward x or y (Vx,y = Vy,x)
          edotxy = pt5 * 
     &            (shdx(iel,1,5)*vl(2,1) + shdx(iel,2,5)*vl(2,2)
     &           + shdx(iel,3,5)*vl(2,3) + shdx(iel,4,5)*vl(2,4)
     &           + shdy(iel,1,5)*vl(1,1) + shdy(iel,2,5)*vl(1,2)
     &           + shdy(iel,3,5)*vl(1,3) + shdy(iel,4,5)*vl(1,4))

c.... second invariant of the strain rate

          edot2 = edotxx * edotxx + edotyy * edotyy
     &            +  two * edotxy * edotxy

          secinv = dsqrt(pt5 * edot2) + 1.0e-9

          eta_star = 1.0e-5
          sigma_y = tcon(2,mat(iel))
          evisc2(intp) = eta_star + sigma_y/secinv

c.... total viscosity (1/diff_creep + 1/disloc_creep)^(-1)

          evisc(intp) = two/(one/evisc1(intp) + one/evisc2(intp)) 
        enddo
c********************************************************************      

      elseif ( ntimvs .eq. 6 ) then
c
c..... Frank Kamenetskii type appoximation
c
         do intp = 1, nipt
            tq = 0.d0
            zq = 0.d0
            do i = 1,4
C     make non-dimensional again
               tq = tq + (tl(i) - Tfac)   * shl(i,intp)
c     depth
               zq = zq + xl(2,i) * shl(i,intp)
            enddo
c
c.... viscosity \eta = visc_0 *exp(E(0.5-T) + V(1-z)
c     where T and z are non-dimensionalized and
c
c     E  = tcon(1,*)
c     V  = tcon(2,*)
c

            evisc(intp) = visc(mat(iel)) *
     &           exp( tcon(1,mat(iel)) * (0.5d0 - tq) +
     &                tcon(2,mat(iel)) * (1.0d0 - zq) )

            if (evisc(intp) .gt. vismax) then
               evisc(intp) = vismax
            endif
            
            if (evisc(intp) .lt. vismin) then
               evisc(intp) = vismin
            endif
c            if((zq.gt.0.45).and.(zq.lt.0.55))then
c               print *,tcon(1,mat(iel)),tq,tcon(2,mat(iel)),zq,
c     &              evisc(intp)
c             endif
         enddo
      endif
c********************************************************************      
      return
c
      end
c
c********************************************************************      
c
      double precision function visexp( t, estar, z, vstar )
c
c This is where the temperature-dependence is set.  
c Notice in this formulation, estar and vstar are 
c dimensional and that a slew of dimensional constants
c are set in this routine.   While this goes against
c the philosophy of ConMan in other scalings, it seems
c to make more sense for viscosity laws.
c sdk 9/19/93
c
c  Input Estar in kJ/mol and Vstar in cm**3/mol
c
      implicit double precision(a-h,o-z)
c
      include 'common.h'
c
c     Kelvin
      Tmin = 273.0d0
      tn = t * 2000.0
c     Joules/(mole*kelvin)
      R    = 8.3144d0
c     kg/m^3
      rho  = 3.3d3
c     m/s^2      
      g    = 9.8
c     1.0d-6 converts cm**3 to m**3
      dscale = rho * g * 1.0d-6
      press = (one - z) * dscale * 2.89d6 
      d0    = dscale * 2.89e6 
c
c This is scaled s.t. at T=1 and depth = 2890 km visexp = 1.0
c 1.0d3 converts kJ/mole to J/mole 
c
c values corresponding to estar and vstar from Karato & Wu (1993) table 1
c wet diffusive olivine (vals next line)
c A 5.3e15 (s^-1) ; n 1.0; m 2.5; estar 240 (kJ/mol); vstar 5 (cm^3/mol);
c dry diffusive olivine (vals next line)
c A 8.7e15 (s^-1) ; n 1.0; m 2.5; estar 300 (kJ/mol); vstar 6 (cm^3/mol);      
c
      visexp = (estar*1.e3 + vstar*press)/(R*(tn+Tmin))
c
      return
      end

