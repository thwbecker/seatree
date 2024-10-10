      subroutine rheol( x , v , t , ivel , evisc)  
c
c
c----------------------------------------------------------------------
c
c  This routine calculates the viscosities which are passed to the
c  stiffness matrix in f_vstf.f through the array local array evisc.
c
c  Code is still in developmental stage. Stress-dependent section has
c    been successfully tested with a simple shear problem. Constant and 
c    Temperature-dependent viscosity sections yield the same results as 
c    previous versions of this subroutine. (SSS 6/93)
c
c  Stress-dependent section assumes that stress exponent n = 3.
c
c  Viscosity scaling for non_Newtonian section:
c       low stress => Newtonian rheology
c       high stress => non-Newtonian rheology
c
c----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
#include<common.h>
#include<materials.h>
#include<fe.h>
c
      dimension  x(nsd,*)     , v(ndof,*)        , t(*)       ,
     &           evisc(5)     , xl(2,4) , vl(8), tl(4), det(5),
     &      shdx(4,5), shdy(4,5), strtmp(6)
c
      real*8 maxsecinv, secinv

       maxsecinv = -1.0d6
c
c..... constant viscosity - Newtonian fluid
c
      if ( ntimvs .eq. 0 ) then
        do intp = 1, nipt
          evisc(intp) = visc(mat(ivel))
        enddo
      endif
c
c..... temperature dependent rheology
c
      if ( ntimvs .eq. 1 ) then
c
        do n = 1, nen
          tl(n) = t(    ien(ivel,n) ) 
          xl(2,n) = x( 2, ien(ivel,n) ) 
        enddo
        do intp = 1, nipt
          tq = shl(1,intp) * tl(1)
     &       + shl(2,intp) * tl(2)
     &       + shl(3,intp) * tl(3)
     &       + shl(4,intp) * tl(4)
          zq = shl(1,intp) * xl(2,1)
     &       + shl(2,intp) * xl(2,2)
     &       + shl(3,intp) * xl(2,3)
     &       + shl(4,intp) * xl(2,4)
          power = visexp(tq, tcon(1,mat(ivel)),
     &                   zq, tcon(2,mat(ivel)) )
          evisc(intp) = 1.0e-6*exp ( power )
c         if (evisc(intp) .gt. 1.0e5) evisc(intp) = 1.0e5
c         if (evisc(intp) .lt. 0.0001) evisc(intp) = 0.0001
          evisc(intp) = visc(mat(ivel))*evisc(intp)
c         write(6,*) "evisc(",intp,")=",evisc(intp) , "temp"
        enddo
c
c..... localize coordinates and velocity
c
        do n=1,nen
          xl(1,n) = x(1,ien(ivel,n) ) * 3.0e6
          xl(2,n) = x(2,ien(ivel,n) ) * 3.0e6
          vl(2*n-1) = v(1,ien(ivel,n) ) *1.0e-6/3.0e6
          vl(2*n  ) = v(2,ien(ivel,n) )  *1.0e-6/3.0e6
        enddo
c     
c.... call the global shape function
c        
        call genshg (shdx, shdy , det , xl , ivel)
c
c..... find velocity derivative at element centers
c
c  normal stress x
        strtmp(1) = 1.0*(shdx(1,5)*vl(1) + shdx(2,5)*vl(3)
     &                 + shdx(3,5)*vl(5) + shdx(4,5)*vl(7))
c  normal stress y
        strtmp(2) = 1.0*(shdy(1,5)*vl(2) + shdy(2,5)*vl(4)
     &                 + shdy(3,5)*vl(6) + shdy(4,5)*vl(8))
c  shear stress
        strtmp(3) = 0.5*(shdx(1,5)*vl(2) + shdx(2,5)*vl(4)
     &                 + shdx(3,5)*vl(6) + shdx(4,5)*vl(8)
     &                 + shdy(1,5)*vl(1) + shdy(2,5)*vl(3)
     &                 + shdy(3,5)*vl(5) + shdy(4,5)*vl(7))
c pressure
        strtmp(4) = ( strtmp(1) + strtmp(2) )

c second invariant
        strtmp(6) = ( strtmp(2)*strtmp(2)  
     &              + strtmp(1)*strtmp(1)  
     &              + strtmp(3)*strtmp(3) )

c
c..... scale by material parameters use last time steps viscosity!
c
        if (strtmp(6) .gt. zero ) then 
           secinv = sqrt(0.5*strtmp(6))
        else
           secinv = 1.0
        endif
        
        if (secinv .gt. maxsecinv) then
           maxsecinv = secinv
        endif
c
        xn=3.0
c get in terms of stress, then divide by shear modulus (Karato and Wu, 1992)
        secinv = secinv*1.0e21/tcon(3,mat(ivel))
c       write(6,*) "secinv =", secinv, secinv**((1-xn)/xn)
        do 2600 n = 1, 5
           evisc(n)=evisc(n)*(secinv**((1-xn)/xn))
2600    continue      
c
c.....save this iteration's viscosity for next iterations stress weighting
c
        do 3200 n = 1, 5
          if (evisc(n) .gt. 1.0e3) evisc(n) = 1.0e3
          if (evisc(n) .lt. 0.1) evisc(n) = 0.1
c         evisc(n) = 1.0/(1.0/evisc(n) + 1.0e-5)
3200    continue      
c
c.....end stress dependent viscosity
c
      end if
c
      if ( FIRST .eq. .true.) then
        do intp = 1, nipt
          evisc(intp) = 1.0
        enddo
      endif
c
      return
      end
c
c****************************************************************
c
      double precision function visexp( t, estar, z, vstar )
c
c This is where the temperature-dependence is set.  
c Notice in this formulation, estar and vstar are 
c dimensional and that a slew of dimensional constants
c are set in this routine.   While this goes against
c the philosophy of ConMan in other scalings, it seems
c to make more sense for viscosity laws.
c
      implicit double precision(a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      include 'common.h'
      real*8 n
c
      n = 3.0
      tn   = 3000.0d0*t
      Tmin = 273.0d0
      R    = 8.3154d0
      rho  = 3.3d3
      g    = 9.8
c 1.0d-6 converts cm**3 to m**3
      dscale = rho * g * 1.0d-6
c
      depth = (1.0 - z) * dscale * 2.89d6 
      d0    = dscale * 2.89e6 
c
c 1.0d3 converts kJ/mole to J/mole 
c
      visexp = (estar*1.e3 + vstar*depth)/(n*R*(tn+Tmin))
c
      return
      end
