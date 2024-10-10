      subroutine input
c
c----------------------------------------------------------------------
c This routine will input all the necessary data, allocates the
c array storage requirements, and sets up the appropriate parameters
c for the processing.
c
c----------------------------------------------------------------------
c
c
      use bc
      use mstres
      use materials
      use solver
      use variables
      use sub_bench

      implicit double precision (a-h,o-z)
c
      include 'common.h'
c
c.... Start the cpu timer
c
      call timer ('input   ')
c
c.... -------------------->  Global Information  <--------------------
c
c.... Input the global information

      read (iin,1000) ititle
c
c.... Append the date to ititle

      call mydate()

1     read (iin,*,err=1,end=999) numnp , nelx  , nelz  , iflow , 
     &                           necho , inrstr, iorstr, nodebn, 
     &                           ntimvs, nwrap , itype , isolve
      write (iout,1100) ititle , numnp , nelx  , nelz  , iflow , 
     &                           necho , inrstr,iorstr , nodebn, 
     &                           ntimvs, nwrap , itype , isolve
c
c ConMan is always 2D
c
      nsd = 2
      ndof = 2

c
c.... Allocate memory and set up the global arrays
c
      ntotal = 0
      allocate (x(nsd,numnp), v(ndof,numnp), t(numnp), tdot(numnp), 
     &          STAT = IERR)
      if (IERR .ne. 0) then
        write(*,*) 'allocation of x, v, t and/or tdot fails'
        STOP      
      endif
      do j = 1, numnp
        t(j) = zero ; tdot(j) = zero
        do i = 1, nsd
          x(i,j) = zero
        enddo
      enddo
      do i = 1, ndof
        do j = 1, numnp
          v(i, j) = zero
        enddo
      enddo
c
c.... if itype = 5, then set wedge to true and use BA
c
      sb_is_wedge_benchmark = 0
      if (itype .eq. 5) then
        itype = 4
        sb_is_wedge_benchmark=1
      endif
c
c     don't want to mess with the 5 switch, though this could probably be combined
c
      if(itype.eq.6)then
         sb_use_tmovis_for_uplate = 1 ! whatever is read in for tmovic will be Batchelor speed
         sb_is_wedge_benchmark=1
         itype=4
      endif
      
      ntotal = ntotal + numnp * (nsd+ndof+2)*2     
c
c.... Read the time integration information

      time = zero
2     read (iin,*,err=2,end=999) nstep ,accel
c
c.... automatically pick the correct parameters based on the temp solver
c
      niter=1
      alpha=1.0
      if (isolve .eq. 1) then
         niter=2
         alpha=0.5
      endif
      write (iout,1200) ititle,  nstep ,niter, alpha, accel
c
c.... Read output information card
c
3     read (iin,*,err=3,end=999) nstprt ,  tmax  , datasv , tsave , 
     &                           tmovis
      write (iout,1400) ititle,  nstprt ,  tmax  , datasv , tsave , 
     &                           tmovis
      if(sb_use_tmovis_for_uplate.eq.1)then
         sb_uplate=tmovis
         tmovis = tsave
      endif
c     
c.... -------------------->  Nodal Information  <--------------------
c
c.... Read the coordinates

      call genfl (x, nsd)

      if (necho .eq. 1) then
        do 200 n = 1, numnp
          if (mod(n,50) .eq. 1) write (iout,1500) ititle,(j,j=1,nsd)
c this is freaken weird...
c         write (iout,*) n, x(1,n), x(2,n)
          write (iout,1600) n, x(1,n), x(2,n)
200     continue
      endif

c
c     check on the box height for the subduction problem 
c
      if(x(2,1).gt.x(2,nelz+1))then
         sb_box_height = x(2,1)
      else
         sb_box_height = x(2,nelz+1)
      end if
c
c.... Read the Boundary Codes in BC array

      allocate (idv(ndof,numnp), idt(numnp), STAT = IERR)
      if (IERR .ne. 0) then
        write(*,*) 'allocation of idv and/or idt fails'
        STOP      
      endif
      do j = 1, numnp
        idt(j) = 0
        idv(1,j) = 0
        idv(2,j) = 0
      enddo 
          
      ntotal = ntotal + numnp*(ndof+1)
      call geniBC (idv,ndof,numnp,neqv,necho)
      ione=1
      call geniBC (idt,ione,numnp,neqt,necho)
c
c.... Read the boundary conditions
c
      allocate (vbc(ndof,numnp), tbc(numnp), STAT = IERR)
      if (IERR .ne. 0) then
        write(*,*) 'allocation of vbc and/or tbc fails'
        STOP      
      endif
      do j = 1, numnp
        tbc(j) = zero
        vbc(1,j) = zero
        vbc(2,j) = zero
      enddo 

      ntotal = ntotal + numnp*(ndof+1)*2
      call genBC (vbc,ndof,numnp,1,necho,x,idv)
      call genBC (tbc,ione,numnp,2,necho,x,idv)
c
c.... Read the boundary nodes for the nusselt number smoother
c
      allocate (nb(2,nodebn), STAT = IERR)
      if (IERR .ne. 0) then
        write(*,*) 'allocation of nb fails'
        STOP      
      endif
      do i = 1, nodebn
        nb(1, i) = 0 ; nb(2, i) = 0
      enddo

      ntotal = ntotal + nodebn*2
c
c.... elements and nodes for Nusselt Number calculation
c
      call genNBC (nb) 
c
c.... read initial temperature 
c
4     read (iin,*,err=4,end=999) pert , xsize , zsize  
      write(iout,1800) ititle, pert, xsize, zsize
c
c.... read equation of state info
c
5     read (iin,*,err=5,end=999) Di, T0, diff_T, cgamma, rho0
      write(iout,1900) ititle, Di, T0, diff_T, cgamma, rho0

c
c
c.... --------------->  Element Group Information  <---------------
c
c.... Read ELEMENT GROUP data

6     read (iin,*,err=6,end=999) numat , numsuf 
c
c.... these used to be read in but never change
c
      numel = nelx*nelz
      nen = 4
      nipt = 5
c
c.... print the element control parameter
c
      write (iout,1700) ititle, numel, numat , numsuf 

c
c.... input element data
c
        call elminp()
c
c.... calculate parent shape element, global shape elements
c     are now calculated on the fly, to save memory.

        call genshg(x)
c
c.... calculate equation of state
c
        call deos(x, t, dens, tref)
c
c.... set the initial temperature and bouyancy field if applicable
c
        call intemp(x, t, tdot)
c
        if (isolve .eq. 1) then              ! Explicit only
c.... from the thermal mass lumped matrix
         call form_temp_mass_matrix( x , tlhs )          
        endif
c
c.... Read surface force / flux boundary conditions
c
      if (numsuf .gt. 0) then
        allocate (isuf(numsuf,2), suf(numsuf,3), STAT = IERR)
        do i = 1, numsuf
          isuf(i,1) = 0    ; isuf(i,2) = 0
          suf(i,1)  = zero ; suf(i,2)  = zero ; suf(i,3) = zero
        enddo
        if (IERR .ne. 0) then
          write(*,*) 'allocation of isuf and/or suf fails'
          STOP      
        endif      
        call rdflux()
      end if          
c
c.... return

      return
c
c.... end of file error handling

999   call error ('input   ','end file',iin)

1000    format(80a1)
1100    format(' ',80a1,//,
     & ' C o n t r o l__ I n f o r m a t i o n                ',  //,
     & ' number_of_nodal_points .....................(numnp_)=',i10//,
     & ' number_of_elements_in_the_x1_direction......(nelx__)=',i10//,
     & ' number_of_elements_in_the_x2_direction......(nelz__)=',i10//,
     & ' execution_code..............................(iflow_)=',i10//,
     & '    eq. 0, data_check_only                            ',  // ,
     & '    eq. 1, execution                                  ',  // ,
     & ' input_echo_parameter........................(necho_)=',i10//,
     & ' input_restart_parameter.....................(inrstr)=',i10//,
     & ' output_restart_parameter....................(iorstr)=',i10//,
     & ' number_of_boundary_nodes_for_nusselt_smooth.(nodebn)=',i10//,
     & ' temperature_dependent_viscosity.............(ntimvs)=',i10//,
     & ' number_of_nodes_for_wrap_around_boundary....(nwrap_)=',i10//,
     & ' form(5=BC+wedge,4=BC,3=EBA,2=TALA,1=ALA)....(itype_)=',i10//,
     & ' temp solver(1=explicit,3=Picard,2=implicit).(isolve)=',i10//)

1200  format(' ',80a1,//,
     & ' T i m e__I n t e g r a t i o n__I n f o r m a t i o n',//,
     & '______Nstep______Niter______Alpha______accel',//,
     &     1x,i10,i8,2x,f10.5,2x,1p,e15.5,//)
1400  format(' ',80a1,//,
     & ' O u t p u t___I n f o r m a t i o n                  ',1p,//,
     & ' number_of_time_steps_per_output_............(nstprt)=',i10//,
     & ' time_of_the_experiment_runs.................(tmax__)=',e14.8//,
     & ' interval_for_flux_and_compressible data.....(datasv)=',e14.8//,
     & ' interval_of_time_output_for_GMT_processing..(tsave_)=',e14.8//,
     & ' interval_for_movie_generation...............(tmovis)=',e14.8//
     &      )

1500  format(' ',80a1,//,
     & ' N o d a l__C o o r d i n a t e s                  ',//,
     & '    Node     ',12x,3('x',i1,:,17x))
1600  format(1p,2x,i5,13x,3(1pe12.5,7x))
1700  format(' ',80a1,//,
     & ' E l e m e n t__G r o u p__C o n t r o l__P a r a m ',//,
     & ' number_of_elements..........................(numel_)=',i10//,
     & ' number_of_materials.........................(numat_)=',i10//,
     & ' number_of_surface_forces....................(numsuf)=',i10//)
1800  format(' ',80a1,//,
     & ' I n i t i a l__T e m p e r a t u r e__P a r a m    ',//,
     & ' perturbation...............................(pert__)=',e14.8//,
     & ' horizontal dimension.......................(xsize_)=',e14.8//,
     & ' vertical dimension.........................(zsize_)=',e14.8//)
1900      format(' ',80a1,//,
     & ' E q u a t i o n__o f__S t a t e__P a r a m e t e r s ',//,
     & ' Dissipation number.........................(Di____)=',e14.8//,
     & ' Surface Temperature (dimensional kelvin)...(T0____)=',e14.8//,
     & ' Temperature difference (dimensional kelvin)(diff_T)=',e14.8//,
     & ' Gruneisen parameter (non-dimensional)......(cgamma)=',e14.8//,
     & ' Reference Density (non-dimensional)........(rho0__)=',e14.8//)

      end
