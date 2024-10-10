        subroutine time_driver ( )
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
      use bc
      use materials
      use fe
      use mstres
      use solver
      use variables

      implicit double precision (a-h,o-z)
c
      include 'common.h'
c
      logical FACSTW
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
c.... density and thermal expansivity

      call deos(x, t, dens, tref) 
      call ddiff(x, t, diff)

c
c.... form the velocity stiffness matrix and residual
c.... STIFFNESS tensor :
c.... T and P dependent viscosity : assemblied per each time step
c.... FACSTW = .true.
c.... Constant viscosity : assemblied only one time
c.... FACSTW = .true. after one step, it becomes .false.
c
c
c if there are no velocities to solve, then don't
c
      if (neqv .ne. 0) then
        call timer('f_vStiff')
        call form_velocity_stiffness_matrix(
     &              x     , v    , t    , idiag,
     &              vlhsa , vlhsb, vbcr , vrhs )   ! asymmetric matrix factorization

        call timer('factor  ')
        call unfact(vlhsa, vlhsb, idiag, neqv)   ! asymmetric matrix solver
c
c add thermal bouyancy to velocity rhs
c
        call timer('f_vRes  ')
        if (itype .eq. 1) then     ! ALA only  
          call dyna_P(x, v, t)
        endif
        call form_velocity_rhs(x, t, vrhs, vbcr)
        call timer('back-slv')
        call unback(vlhsa, vlhsb, vrhs, idiag, neqv)  ! asymmetric matrix solver
        call vadd  (vrhs , v  )
      endif
c
      if (necho .eq. 1) then
        call print(x,v,t,istep0)
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

        if (time .le. tmax) then
c
c.... pick the correct time step
c
        call  timer('pick_dt ')
        call  pickdt( x , v )

        dt = dt * accel
        time = time + dt

c.... ---------------> predictor phase <---------------

        call timer('predict ')
        fac = (one-alpha)*dt
        do 1000 i = 1 , numnp
          t(i) = t(i) + fac * tdot(i)
          tdot(i) = zero
1000    continue

c.... ---------------> multi-corrector phase <---------------

        do 2000 npass = 1, niter
c
c.... form the  residual matrix for and solve for temperature
c
          call timer('f_tRes  ')
          call deos(x, t, dens, tref) 
c
c.... first call for temperature equation
c
          if ( (isolve .eq. 2) .or. (isolve .eq. 3) ) then    ! Implicit & Picard
             call form_temp_rhs_implicit(t    , tdot ,
     &                x    , v    , idiagt, trhs , tlhsa, tlhsb)
          elseif (isolve .eq. 1) then                         ! Explicit
             call form_temp_rhs_explicit(t    , tdot ,
     &                x    , v    ,  tlhs , trhs   )
          endif
c
c.... end of the multi-corrector loop
c
2000    continue  

c.... form the mass and residual matrix for velocity

        if (ntimvs .gt. 0) then
          FACSTW  = .true.
        else
          FACSTW  = .false.
        end if
c
c skip velocity solve if there are no unknowns
c
        if (neqv .ne. 0) then
          if (itype .eq. 1) then     ! ALA only  
            call dyna_P(x, v, t)
          endif
          if ( FACSTW ) then               ! Stiffness
            call timer('f_vStiff')
            call form_velocity_stiffness_matrix(
     &             x     , v    , t    , idiag ,
     &             vlhsa , vlhsb, vbcr , vrhs )  

            call timer('factor  ')
            call unfact(vlhsa, vlhsb, idiag, neqv)   ! asymmetric matrix solver
          endif                            ! Stiffness
c
c add thermal bouyancy to velocity rhs
c
          call timer('f_vRes  ')
          call form_velocity_rhs(x, t, vrhs, vbcr)
          call timer('back-slv')
          call unback(vlhsa, vlhsb, vrhs, idiag, neqv)  ! asymmetric matrix solver
          call vadd  (vrhs , v  )
        endif
c
c put this after to see if it would work for the second iteration
c
c output is all in subroutine output () so that user can modify
c
        call timer('output  ')
        lstep = lstep + 1
        call output ()
c
c.... end of the NSTEP loop
c     
        endif
c
3000  continue

      return
      end
