      subroutine proces
c
c----------------------------------------------------------------------
c
c  This subroutine set up the parameters for the global processing
c procedure, and calls up the right global time integration algorithm
c
c----------------------------------------------------------------------
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
#include<common.h>
#include<variables.h>
#include<bc.h>
c
c.... check for execution mode
c
      if (iflow .eq. 0) return
c
c.... engage the timer
c
      call timer ('process ')
c
c.... call the iterative solver--time integration driver routine
c
      call timdrv (x   , v   , t   , tdot  ,
     &             idv , idt , vbc , tbc   )
c
c.... return
c
      return
      end
