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
     &           evisc(5) 
      dimension xl(4) , tl(4)
c
c..... constant viscosity - Newtonian fluid
c
      if ( ntimvs .eq. 0 ) then
        do intp = 1, nipt
          evisc(intp) = visc(mat(ivel))
        enddo
        return
      endif
c
c..... temperature dependent rheology
c
      if ( ntimvs .eq. 1 ) then
c
        do n = 1, nen
          tl(n) = t(    ien(ivel,n) ) 
          xl(n) = x( 2, ien(ivel,n) ) 
        enddo
        do intp = 1, nipt
          tq = shl(1,intp) * tl(1)
     &       + shl(2,intp) * tl(2)
     &       + shl(3,intp) * tl(3)
     &       + shl(4,intp) * tl(4)
          zq = shl(1,intp) * xl(1)
     &       + shl(2,intp) * xl(2)
     &       + shl(3,intp) * xl(3)
     &       + shl(4,intp) * xl(4)
          power = visexp(tq, tcon(1,mat(ivel)),
     &                   zq, tcon(2,mat(ivel)) )
c         evisc(intp) = 1.0e-9*exp ( power ) 
          evisc(intp) = exp ( power ) 
          if (evisc(intp) .gt. tcon(3,mat(ivel))) 
     &           evisc(intp) = tcon(3,mat(ivel))
          if (evisc(intp) .lt. 1.0e-1) evisc(intp) = 1.0e-1
          evisc(intp) = visc(mat(ivel))*evisc(intp)
        enddo
        return
      end if
c
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
c sdk 9/19/93
c
c  Input Estar in kJ/mol and Vstar in cm**3/mol
c  and this is scaled s.t. at T=1 and depth = 150 km visexp = 1.0
c
      implicit double precision(a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      include 'common.h'
c
      tn   = 1100.0d0*t
      Tmin = 773.0d0
      R    = 8.3154d0
      rho  = 3.3d3
      g    = 9.8
c 1.0d-6 converts cm**3 to m**3
      dscale = rho * g * 1.0d-6
c
c     depth = (600.0 - z) * dscale * 2.89d6 
      depth = (1.0 - z) * dscale * 2.89d6 
      d0    = dscale * 2.89e6 
c
c 1.0d3 converts kJ/mole to J/mole 
c
      visexp = (estar*1.e3 + vstar*depth)/(R*(tn+Tmin))
     &       - (estar*1.e3 + vstar*depth)/(R*(1100.0+Tmin))
c
c     if (visexp .gt. 10.0) visexp = 10.0
      return
      end
c
c****************************************************************
c
