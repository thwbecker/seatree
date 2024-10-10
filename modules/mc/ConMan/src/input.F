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
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
#include<common.h>
#include<stress.h>
#include<variables.h>
#include<bc.h>
c
c.... Start the cpu timer
c
      call timer ('input   ')
c
c.... -------------------->  Global Information  <--------------------
c
c.... Input the global information
c
      read (iin,1000) ititle
c
c.... Append the date to ititle
c
      call mydate()
c
1     read (iin,*,err=1,end=999) numnp , nsd   , ndof , nelx  ,nelz  ,
     &                           mprec , iflow , necho, inrstr,iorstr,
     &                           nodebn, ntimvs, ntseq, numeg ,
     &                           isky  , nwrap 
      write (iout,1100) ititle , numnp ,  nsd  , ndof , nelx  ,nelz  ,
     &                           mprec , iflow , necho, inrstr,iorstr,
     &                           nodebn, ntimvs, ntseq, numeg ,
     &                           isky  , nwrap 
c
c.... Allocate memory and set up the global arrays
c
      ntotal = 0
      call mmgetblk('x     ','mm', ipx   ,numnp*nsd , 2 , icode)
      call mmgetblk('v     ','mm', ipv   ,numnp*ndof, 2 , icode)
      call mmgetblk('t     ','mm', ipt   ,  numnp   , 2 , icode)
      call mmgetblk('tdot  ','mm', iptdot,  numnp   , 2 , icode)
      ntotal = ntotal + numnp*(nsd+ndof+2)*2
c
c.... Read the time integration information
c
      time = zero
      if (ntseq .ne. 1) call error('input   ','ntseq   ', ntseq)
2     read (iin,*,err=2,end=999) nstep ,niter, alpha, accel, epstol
      write (iout,1200) ititle,  nstep ,niter, alpha, accel, epstol
c
c.... Read output information card
c
3      read (iin,*,err=3,end=999) nsdprt , nsvprt , nstprt , nsmprt 
       write (iout,1400) ititle,  nsdprt , nsvprt , nstprt , nsmprt 
c
c.... -------------------->  Nodal Information  <--------------------
c
c.... Read the coordinates
c
      call genfl (x, nsd)
c
      if (necho .eq. 1) then
        do 200 n = 1, numnp
          if (mod(n,50) .eq. 1) write (iout,1500) ititle,(j,j=1,nsd)
c this is freaken weird...
          write (iout,*) n, x(1,n), x(2,n)
c         write (iout,1600) n, x(1,n), x(2,n)
200     continue
      endif
c
c.... Read the Boundary Codes in BC array
c
      call mmgetblk('idv   ','mm', ipidv ,numnp*ndof, 1 , icode)
      call mmgetblk('idt   ','mm', ipidt , numnp    , 1 , icode)
      ntotal = ntotal + numnp*(ndof+1)
      call geniBC (idv,ndof,numnp,neqv,necho)
      call geniBC (idt,1   ,numnp,neqt,necho)
c
c.... Read the boundary conditions
c
      call mmgetblk('vbc   ','mm', ipvbc ,numnp*ndof, 2 , icode)
      call mmgetblk('tbc   ','mm', iptbc , numnp    , 2 , icode)
      ntotal = ntotal + numnp*(ndof+1)*2
      call genBC (vbc,ndof,numnp,1,necho,zero,x,idv)
      call genBC (tbc, 1  ,numnp,2,necho,zero,x,idv)
c     call genBC (vbc,ndof,numnp,1,necho,zero)
c     call genBC (tbc, 1  ,numnp,2,necho,zero)
c
c.... Read the boundary nodes for the nusselt number smoother
c
      call mmgetblk('nb    ','mm', ipnb  , nodebn*2 , 1 , icode)
      ntotal = ntotal + nodebn*2
      call genNBC (nb) 
c
c.... set the initial temperature and bouyancy field if applicable
c
      call intemp(x, t, tdot)
c
c
c.... --------------->  Element Group Information  <---------------
c
c.... Read ELEMENT GROUP data
c
4     read (iin,*,err=4,end=999) ntype , numel ,    nen , nenl ,
     &                           numat , nedof , numsuf , nipt
c
       nEGdf = ndof * numnp
c temporary definition
       nEGnp = numnp
c
c.... print the element control parameter
c
      write (iout,1700) ititle, ntype, numel, numat , numsuf , nen, 
     &                  nipt  , nEGnp, nEGdf
c
      call EGlib ('input   ')
c
c.... Read surface force / flux boundary conditions
c
      if (numsuf .gt. 0) then
        call mmgetblk('isuf  ','mm', ipisuf, numsuf*2 , 1,icode)
        call mmgetblk('suf   ','mm', ipsuf , numsuf*3 , 2,icode)
        call rdflux()
      end if          
c
c.... return
c
      return
c
c.... end of file error handling
c
999   call error ('input   ','end file',iin)
c
1000    format(80a1)
1100    format(' ',80a1,//,
     & ' C o n t r o l   I n f o r m a t i o n                ',  //,
     & ' number of nodal points  . . . . . . . . . . (numnp )=',i10//,
     & ' number of space dimensions  . . . . . . . . (nsd   )=',i10//,
     & ' number of degrees of freedom per node . . . (ndof  )=',i10//,
     & ' number of elements in the x1 direction. . . (nelx  )=',i10//,
     & ' number of elements in the x2 direction. . . (nelz  )=',i10//,
     & ' single or double precision  . . . . . . . . (mprec )=',i10//,
     & ' execution code  . . . . . . . . . . . . . . (iflow )=',i10//,
     & '    eq. 0, data check only                            ',  // ,
     & '    eq. 1, execution                                  ',  // ,
     & ' input echo parameter  . . . . . . . . . . . (necho )=',i10//,
     & ' input restart parameter   . . . . . . . . . (inrstr)=',i10//,
     & ' output restart parameter  . . . . . . . . . (iorstr)=',i10//,
     & ' number of boundary nodes for nusselt smooth (nodebn)=',i10//,
     & ' temperature dependent viscosity . . . . . . (ntimvs)=',i10//,
     & ' number of time sequences  . . . . . . . . . (ntseq )=',i10//,
     & ' number of element groups  . . . . . . . . . (numeg )=',i10//,
     & ' flag for skyline or frontal solve . . . . . (isky  )=',i10//,
     & ' number of nodes for wrap around boundary. . (nwrap )=',i10//)

1200  format(' ',80a1,//,
     & ' T i m e   I n t e g r a t i o n   I n f o r m a t i o n',//,
     & '       Nstep    Niter    Alpha     accel         Eps-tol',//,
     &     1x,i10,i8,2x,f10.5,2x,1p,e15.5,1x,e15.5,//)
1400  format(' ',80a1,//,
     & ' O u t p u t   I n f o r m a t i o n                  ',1p,//,
     & ' number of time steps per output disk. . . . (nsdout)=',i10//,
     & ' number of time steps per output velocity. . (nsvprt)=',i10//,
     & ' number of time steps per output temperature (nstprt)=',i10//,
     & ' number of time steps per output mean. . . . (nsmprt)=',i10//)
1500  format(' ',80a1,//,
     & ' N o d a l   C o o r d i n a t e s                  ',//,
     & '    Node     ',12x,3('x',i1,:,17x))
1600  format(1p,2x,i5,13x,3(1pe12.5,7x))
1700      format(' ',80a1,//,
     &  ' E l e m e n t   G r o u p   C o n t r o l   P a r a m ',//,
     &  ' element type  . . . . . . . . . . . . . . . (ntype )=',i10//,
     &  '    eq. 2, 2D element                                 ',  / ,
     &  '    eq. 3, 3D element                                 ',  //,
     &  ' number of elements  . . . . . . . . . . . . (numel )=',i10//,
     &  ' number of materials . . . . . . . . . . . . (numat )=',i10//,
     &  ' number of surface forces  . . . . . . . . . (numsuf)=',i10//,
     &  ' number of element nodes . . . . . . . . . . (nen   )=',i10//,
     &  ' number of integration points  . . . . . . . (nipts )=',i10//,
     &  ' number of group nodal points  . . . . . . . (nEGnp )=',i10//,
     &  ' number of group dofS  . . . . . . . . . . . (nEGdf )=',i10//)
      end
