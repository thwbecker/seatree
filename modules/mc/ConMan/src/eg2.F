        subroutine EG2 (task)
c
c----------------------------------------------------------------------
c This program is the Element Group 2D routine
c
c input:
c  npar (*)      : the parameters of the current element group
c  npt  (*)      : pointers to the element group arrays
c  task          : task to be performed
c
c----------------------------------------------------------------------
c
c
      implicit double precision (a-h,o-z)
      dimension psol(111,101)
c
c.... deactivate above card(s) for single precision operation
c
#include<common.h>
#include<solver.h>
#include<stress.h>
#include<variables.h>
#include<materials.h>
c
      character*8 task
      logical FACSTW
c
c.... ---------------------->  'input   '  <----------------------
c
      if (task .eq. 'input   ') then
c
c.... input element data 
c
        call elminp()
c
c.... calculate parent shape element, global shape elements
c      are now calculated on the fly, to save memory.
c
        call genshp ( x )
#ifdef IMPLICIT
#else
c
c.... from the thermal mass lumped matrix
c
        call tmass ( x , tlhs ) 
#endif
        return
      endif
c
c.... ---------------------->  'pick_dt '  <----------------------
c
c.... pick the timestep based on the courant time step
c
        if (task .eq. 'pick_dt ') then
          call  timer('pick_dt ')
          call  pickdt( x , v )
          return
        end if
c
c.... ---------------------->  'v_mat_Rs'  <----------------------
c
c.... form the velocity stiffness matrix and residual
c
      if (task .eq. 'v_mat_Rs') then
c
        if ( FACSTW ) then
          call timer('f_vStiff')
          call f_vStf(x    , v    , t    ,  idiag    , 
     &                vlhs , vbcr , vrhs )
c
          call timer('factor  ')
          if (isky .eq. 1) then
            call factor (vlhs , idiag   , neqv,  vrhs)
          elseif (isky .eq. 0) then
            call vfactor(vlhs , idiag   , neqv,  vrhs)
          else
            call error('factor  ','isky =  ',isky  )
          endif
c 
        endif
        call timer('f_vRes  ')
c
c add thermal bouyancy to velocity rhs
c
        call f_vRes(x , t , vrhs , vbcr)
c
        call timer('back-slv')
        if (isky .eq. 1) then
          call back  (vlhs , vrhs ,  idiag    , neqv) 
        elseif (isky .eq. 0) then
          call vback (vlhs , vrhs ,  idiag    , neqv) 
        else
          call error('back-slv','isky =  ',isky  )
        endif
        call vadd   (vrhs , v )
c
        return
      endif
c
c.... ---------------------->  'RES_temp'  <----------------------
c
c.... form the temperature  residual and solve for current temperature
c
      if (task .eq. 'RES_temp') then
c
        call timer('f_tRes  ')
c
c.... first call for temperature equation
c
#ifdef IMPLICIT
        call f_tRes(t    , tdot ,
     &              x    , v    , idiagt ,
     &              trhs , tlhsa, tlhsb     )
#else
        call f_tRes(t    , tdot , 
     &              x    , v    , idiag  ,
     &              tlhs , trhs     ) 
#endif
        return
      endif
c
c.... ----------------------> 'prt_flux' <---------------------------
c
      if (task .eq. 'prt_flux') then
c
c trhs is being used as a temp variable for the nusselt number at the node
c
        call fluxke(x   , v   , t  , trhs ) 
        return
      endif
c
c.... ----------------------> 'prt_strn' <---------------------------
c
      if (task .eq. 'prt_str ') then
        call stres ( x , v , t ) 
        call prtstr( x , stress )
        call geoid ( x, t, stress)
c
c trhs is being used as a temp variable for the pressure at the node
c
c         do n=1,numnp
c           dum(n) = stress(4,n)
c         enddo
c         call projct( dum , psol )
c         iflag = 1
c         call prtstr_reg( psol , iflag )
c second invariant
c         do n=1,numnp
c           dum(n) = log10(stress(5,n))
c         enddo
c         call projct( dum , psol )
c         iflag = 0
c         call prtstr_reg( psol , iflag )
c viscosity
c         do n=1,numnp
c           dum(n) = log10(stress(6,n)/(1.0e0))
c         enddo
c         call projct( dum , psol )
c         iflag = 0
c         call prtstr_reg( psol , iflag )
        return
      endif
c
c.... error
c
      call error ('EG2     ', task, 0)
      return
c
      end