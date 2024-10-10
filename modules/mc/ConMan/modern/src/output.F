        subroutine output ( )
c----------------------------------------------------------------------
c
c  All output moved here for convenience 
c
c working arrays:
c  x    (ndof,numnp)         : Global coordinates
c  v    (ndof,numnp)         : Velocity variables 
c  t    (numnp)              : Temperature variables 
c  tdot (numnp)              : Time derivative of temperature variables
c  idv  (ndof,numnp)         : Equation numbers for velocity
c  idt  (numnp)              : Equation numbers for temperature
c  vbc  (ndof,numnp)         : Boundary values for velocity
c  tbc  (numnp)              : Boundary values for temperature
c
c----------------------------------------------------------------------
c
c
      use bc
      use materials
      use mstres
      use solver
      use variables

      implicit double precision (a-h,o-z)
      parameter (ipsolx = 301, ipsolz = 101)
c
      include 'common.h'
c
      logical FACSTW, restart
      dimension t_adia(numnp)
c     dimension tpara1(3,2), npara2(2)
      dimension psol(ipsolx,ipsolz), ipsol(ipsolx,ipsolz)

      character*8 vfile

c
c.... initialization
c
      do i = 1, ipsolx
        do j = 1, ipsolz
          psol(i,j) = zero
          ipsol(i,j) = 0
        enddo
      enddo
c
c.... print out the boundary flux (Nu #) & kinetic energy (in fluxke.F)
c.... print out data corresponding compressible flow such as <phi> and <W>
c.... according to J&M (1980) (in print_compbm_data.F)
c
c     print out data per datasv
c
       if (isolve .ne. 1) then      ! always output for implicit or picard solvers
          call fluxke(x   , v   , t  , trhs )
          call print_compbm_data(x,v,t,lstep)
       elseif ((mod(time, datasv).lt.dt).or.(lstep.eq.nstep)) then ! every datasv for explicit
          call fluxke(x   , v   , t  , trhs )
          call print_compbm_data(x,v,t,lstep)
       endif         
c
       restart = .false.
!=============================================================================== !
c
c.... print out temperature and velocity fields for GMT processing
c
        if ((mod(time, tsave) .le. dt)
     &      .or.(lstep .eq. nstep)
     &      .or.(mod(lstep,nstprt) .eq. 0)) then

c.... print out x(numnp,*), v(numnp,*), t(numnp), t_addi(numnp)
          call print(x,v,t,lstep,time)
          call stres(x,v,t)
          call prtstr(x,stress)
c         call output_rheol(x,v,t)
c      
c.... restart output file every nsdout timesteps      
c
          write(irsout,1001) lstep,time
          write(irsout,1002) 
          do n=1,numnp
            write(irsout,1003) n,t(n),tdot(n)
          enddo
        endif
c
c.... print out movies 
c
        if (mod(time,tmovis) .lt. dt) then
          vfile = 'mov'  ! head of movie files
c.... add adiabatic temperature profile to the potential T
c.... the temperature used in the movie is Celcius is calculated by
c.... removing the surface temperature
          do i = 1, numnp
            t_adia(i) = t(i) +  tref(i) - T0/diff_T
            if (t_adia(i) .gt. 1.0d0) t_adia(i) = 1.0d0
          enddo
c         call grdpar(x, tpara1, npara2)
c         call projct(t_adia, psol, ipsol)
c         call ppmprt(ipsol, npara2, vfile, nmplot)
c         nmplot = nmplot + 1
c
c print mass flux out every movie time
c
c         call masflx()
c
        endif           
c
      return
1001  format('step = ',i7,'___time = ',f10.6)
1002  format('_node_____temperature______tdot' )            
1003  format(1x,i6,4(1pe15.8,1x))
      
      end
