      subroutine deos(x, t, dens, tref)
c----------------------------------------------------------------------
c
c  This routine calculates equation of state for density variation
c  
c  input
c  x : node coodinate    (2,numnp) for 2D
c  t : node temperature  (numnp) for 2D
c  
c  output
c  dens : node density    (numnp) for 2D
c  tref : adiabatic temperature (numnp) for 2D (non-dimensionalized with Kelvin)
c----------------------------------------------------------------------
c
      use fe, only : ien
      implicit double precision (a-h,o-z)
c
#include<common.h>
c
      dimension x(nsd,*), t(*), dens(*), tref(*)
      dimension xl(2,4), tl(4), tlr(4), tlw(4)

      if  ( (itype .eq. 1) .or. (itype .eq. 2) ) then    ! ALA and TALA     
c
c.... loop over the elements
c
        do ivel = 1 , numel
c
c.... set up the parameters
c
          do n = 1, nen
            tl(n)       = t(  ien(ivel,n) ) 
            tlr(n)      = tref(ien(ivel,n) )
            tlw(n)      = tl(n) + tlr(n)
            xl(1,n)     = x(1, ien(ivel,n) ) 
            xl(2,n)     = x(2, ien(ivel,n) ) 
          enddo 

c  Adams-Williamson EOS - see Jarvis & Mckenzie, 1980 eq. 15
c  Di = d / Ht                      (Dissipation Number)
c  cgamma = alpha * Ks / (rho * Cp) (Gruneisen Parameter)
c  rho0 = 3700 (kg/m3) : upper mantle density : 3700 (kg/m3)
c  rho is dimensionless value.

          do n = 1, nen
            dens(ien(ivel,n)) = rho0               ! reference density
     &        * exp( Di / cgamma * (one-xl(2,n)))
            tref(ien(ivel,n)) = T0 / diff_T        ! adiabatic temperature
     &          * (exp( Di * (one -xl(2,n) )))
          enddo
        enddo
      else if  ( (itype .eq. 3) .or. (itype .eq. 4) ) then  ! EBA and BA
        do ivel = 1, numel
          do n = 1, nen
           dens(ien(ivel,n)) = rho0              ! reference density
           tref(ien(ivel,n)) = T0 / diff_T       ! adiabatic temperature (surface T in K)
          enddo
        enddo
      endif
c
c....      
c
      return
      end
