      subroutine elminp
c
c----------------------------------------------------------------------
c
c  This routine inputs the element data for EG2
c
c
c----------------------------------------------------------------------
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
#include<common.h>
#include<solver.h>
#include<materials.h>
#include<variables.h>
#include<fe.h>
#include<bc.h>
#include<stress.h>
c
c.... read the material properties
c
      call mmgetblk('mat-no ','mm', ipmat , numel   , 1 , icode)
      call mmgetblk('visc   ','mm', ipvisc, numat   , 2 , icode)
      call mmgetblk('alam   ','mm', ipalam, numat   , 2 , icode)
      call mmgetblk('diff   ','mm', ipdiff, numat   , 2 , icode)
      call mmgetblk('ra     ','mm', ipra  , numat   , 2 , icode)
      call mmgetblk('mhu    ','mm', ipmhu , numat   , 2 , icode)
      call mmgetblk('tcon   ','mm', iptcon, numat*3 , 2 , icode)
      ntotal = ntotal + numat*8*2 + numel
      call genmat ()
c
c.... read and set up IEN array
c
      call mmgetblk('ien    ','mm', ipien , numel*nen, 1 , icode)
      call mmgetblk('lmv    ','mm', iplmv,numel*nen*ndof,1,icode)
      call mmgetblk('lmt    ','mm', iplmt , numel*nen, 1 , icode)
      ntotal = ntotal + numel*nen*(ndof+2)  
c
c diag arrays are for the solver
c
      call mmgetblk('idiag  ','mm', ipidiag, neqv   , 1 , icode)
      ntotal = ntotal + neqv
#ifdef IMPLICIT
      call mmgetblk('idiagt ','mm', ipidiat, neqt   , 1 , icode)
      ntotal = ntotal + neqt
      call genien (ien, idiag , mat, lmv , lmt , idv, 
     &             idiagt  , idt )
#else
      call genien (ien, idiag , mat, lmv , lmt , idv)
#endif
c
c.... generate the  element shape functions and 
c.... calculate min element dimension per node
c
      call mmgetblk('shdx   ','mm', ipshdx, nen*nipt , 2 , icode)
      call mmgetblk('shdy   ','mm', ipshdy, nen*nipt , 2 , icode)
      call mmgetblk('det    ','mm', ipdet , nen*nipt , 2 , icode)
      call mmgetblk('shl    ','mm', ipshl , nen*nipt , 2 , icode)
      ntotal = ntotal + nen*nipt*4*2 
c
c solver memory
c
      call mmgetblk('V RHS  ','mm', ipvrhs , neqv   , 2 ,icode)
      call mmgetblk('vbcr   ','mm', ipvbcr , neqv   , 2 ,icode)
      call mmgetblk('T rhs  ','mm', iptrhs , nEGnp  , 2 ,icode)
      ntotal = ntotal + 2*neqv*2 + nEGnp*2
#ifdef IMPLICIT
      call mmgetblk('T lhs a','mm', iptlhsa, nsizet , 2 ,icode)
      call mmgetblk('T lhs b','mm', iptlhsb, nsizet , 2 ,icode)
      ntotal = ntotal + 2*nsizet*2
#else
      call mmgetblk('T lhs e','mm', iptlhs , nEGnp  , 2 ,icode)
      ntotal = ntotal + nEGnp*2
#endif
      call mmgetblk('pmass  ', 'mm', ippmas , numnp , 2 ,icode)
      call mmgetblk('dummy  ', 'mm', ipdum  , numnp , 2 ,icode)
      call mmgetblk('stress ', 'mm', ipstrs ,numnp*6, 2 ,icode)
      call mmgetblk('V lhs  ','mm', ipvlhs , nsize  , 2 ,icode)
      ntotal = ntotal + nsize*2 + numnp*8*2
      call mmprint(6)
      write(6,*) "Memory (Mb) =", ntotal*4e-6
c
c stress and nusselt number arrays
c
      call mmverify()
c
c.... return
c
      return
c
1000      format(' ',80a1,//,
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
1001      format(' ', //,
     &  ' M a t r i x   O p e r a t i o n s   P a r a m e t e r ',//,
     &  ' number of velocity equations  . . . . . . . (neqv  )=',i10//,
     &  ' Matrix Average bandwidth  . . . . . . . . . (nband )=',i10//)
      end
