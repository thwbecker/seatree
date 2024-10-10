        subroutine timdrv (x     , v     , t      , tdot   ,
     &                     idv   , idt   , vbc    , tbc    )
c----------------------------------------------------------------------
c
c  This is the driver routine for Multi-Element Group Predictor 
c Multi-Corrector Algorithm, currently using LDU solver.
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
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      include 'common.h'
c
        logical FACSTW
        dimension x(ndof,1)   , v(ndof,1)   , t(1)       , tdot(1)  , 
     &            idv(ndof,1) , idt(1)      , vbc(ndof,1), tbc(1)   
        dimension psol(111,101)
c.... initialize
c
c.... load the velocity boundary values into the velocity array
c
      do 100 n = 1 , numnp
        v(1,n) = vbc(1,n)
        v(2,n) = vbc(2,n)
100   continue
c
c.... form the velocity stiffness matrix and residual
c     and solve for the initial velocity 
c **** note for now all equation solving in done on the element level****
c
c.... set factor switch to yes
c
      FACSTW  = .true.
      call EGlib( 'v_mat_Rs')
      if (necho .eq. 1) then
        call print(x,v,t,ndof,nelx,nelz,numnp,istep0,time0)
      end if
c
c.... loop through the time sequences
c
      time=time0
      lstep=istep0
c
c.... loop through the time steps
c
      do 3000 istep = 1, nstep
c
c.... pick the correct time step
c
        call EGlib('pick_dt ')

        dt = dt*accel
        time=time+dt
c
c.... ---------------> predictor phase <---------------
c
      call timer('predict ')
      fac = (one-alpha)*dt
      do 1000 i = 1 , numnp
        t(i) = t(i) + fac * tdot(i)
        tdot(i) = zero
1000  continue
c
c.... ---------------> multi-corrector phase <---------------
c
      do 2000 npass = 1, niter
c
c.... form the  residual matrix for and solve for temperature
c

      call EGlib ('RES_temp')
c
c.... end of the multi-corrector loop
c
2000  continue
c
c.... form the mass and residual matrix for velocity
c
      if (ntimvs .eq. 1) then
         FACSTW  = .true.
      else
         FACSTW  = .false.
      end if
c
      call EGlib ('v_mat_Rs')
c
c.... print out the boundary fluxes & kinetic energy
c
      call timer('output  ')
      lstep = lstep + 1
      call EGlib ('prt_flux')
      call mmverify()
c
c.... print out temperature and velocity fields
c
      if (mod(lstep, nstprt) .eq. 0) then
        call print(x,v,t,ndof,nelx,nelz,numnp,lstep,time)

      endif
c
c.... print out stress fields
c
      if (mod(lstep, nsmprt) .eq. 0)then
         call EGlib ('prt_str ')
      endif
c
c.... overwrite restart output file every nsdout timesteps      
c
      if ( iorstr .eq. 1) then
        if (mod(lstep, nsdprt) .eq. 0) then
           open(unit=irsout, file=name5, status='unknown') 
           rewind(irsout)
           write(irsout,1001) lstep,time
           write(irsout,1002) 
           do 9000 n=1,numnp
              write(irsout,1003) n,t(n),tdot(n)
9000       continue
           close(irsout) 
        end if
      end if
c
c.... end of the NSTEP loop
c
3000  continue
c
c.... return
c
1001  format('step = ',i5,'   time = ',f10.6)
1002  format(' node     temperature      tdot' )            
1003  format(1x,i5,4(1pe15.8,1x))
      return
      end
