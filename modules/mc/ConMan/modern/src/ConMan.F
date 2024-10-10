c----------------------------------------------------------------------
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c **                                                                 **
c **              ***** ==>    C O N M A N    <== *****              **
c **                                                                 **
c **                    Convective Mantle Flow code                  **
c **                                                                 **
c **       A Petrov Galerkin Finite Element Program to Analyze       **
c **    incompressible and compressible Thermal Convective Flows.    **
c **                                                                 **
c **                                                                 **
c **   written by:  Arthur Raefsky, Scott King, Brad Hager,          **
c **                       and Changyeol Lee                         **
c **   with contributions by many others documented in the manual.   **
c **                                                                 **
c **                                                                 **
c **       Compatible with GNU fortran 90 or higher version          **
c *********************************************************************
c *********************************************************************
c----------------------------------------------------------------------

c Declaration of global variables
      module bc
        implicit double precision (a-h,o-z)
        pointer :: idv(:,:)     ! idv(ndof, numnp)
        pointer :: vbc(:,:)     ! vbc(ndof, numnp)
        pointer :: idt(:)       ! idt(numnp)
        pointer :: tbc(:)       ! tbc(numnp)
        pointer :: isuf(:,:)    ! isuf(numsuf,2)
        pointer :: suf(:,:)     ! suf(numsuf,3)    
      end module bc
      
      module fe
       implicit double precision (a-h,o-z)
       pointer :: ien(:,:)     ! ien(numel, nen)
       pointer :: lmv(:,:)     ! lmv(numel, nen*ndof)
       pointer :: lmt(:,:)     ! lmt(numel, nen)
       pointer :: shl(:,:)     ! shl(nen, nipt)
       pointer :: shdx(:,:,:)  ! shdx(numel,nen,nipt)
       pointer :: shdy(:,:,:)  ! shdy(numel,nen,nipt)
       pointer :: det(:,:)     ! det(numel,nipt)
        pointer :: work(:)     ! work(nwork)    
      end module fe
   
      module materials
       implicit double precision (a-h,o-z)
       pointer :: dens(:)      ! dens(numnp)  : density
       pointer :: diff(:)      ! diff(numnp)  : thermal diffusivity
       pointer :: tref(:)      ! tref(numnp)  : adiabatic temperature 
       pointer :: mat(:)       ! mat(numel)
       pointer :: visc(:)      ! visc(numat)
       pointer :: ra(:)        ! ra(numat)
       pointer :: alam(:)      ! alam(numat)
       pointer :: dmhu(:)      ! dmhu(numat)
       pointer :: tcon(:,:)    ! tcon(2,numat)
      end module materials

      module solver
       implicit double precision (a-h,o-z)
       pointer :: vlhsa(:)   ! vlhsa(nsize)  upper triangle matrix
       pointer :: vlhsb(:)   ! vlhsb(nsize)  lower triangle matrix
       pointer :: vrhs(:)    ! vrhs(neqv)
       pointer :: vbcr(:)    ! vbcr(neqv)
       pointer :: idiag(:)   ! idiag(neqv)
       pointer :: trhs(:)    ! trhs(numnp)
       pointer :: tlhsa(:)   ! tlhsa(nsizet)  Picard & Imiplicit
       pointer :: tlhsb(:)   ! tlhsb(nsizet)  Picard & Imiplicit
       pointer :: idiagt(:)  ! idiagt(neqt)   Picard & Imiplicit
       pointer :: tlhs(:)    ! tlhs(numnp)    Explicit
      end module solver

      module mstres
       implicit double precision (a-h,o-z)
       pointer :: nb(:,:)    ! nb(2,nodebn)
       pointer :: stress(:,:)! stress(5,numnp)
       pointer :: pmass(:)   ! pmass(numnp)
      end module mstres

      module variables
       implicit double precision (a-h,o-z)
       pointer :: x(:,:)     ! x(nsd,numnp)
       pointer :: v(:,:)     ! v(ndof,numnp)
       pointer :: t(:)       ! t(numnp)
       pointer :: tdot(:)    ! tdot(numnp)
      end module variables

C     subduction benchmark parameters
      module sub_bench
       integer :: sb_use_tmovis_for_uplate = 0 ! use tmovis to replace sb_uplate
       integer :: sb_is_wedge_benchmark = 0 ! should we check if in fault etc.
       double precision :: sb_box_height = 600.d0 ! default value for
                                                  ! box heigh, will get
                                                  ! adjusted
       double precision :: sb_plate_thickness = 50.d0 ! stagnant lid
                                                      ! thickness
       double precision :: sb_uplate= 2.178773914605951d0  ! 5
                                                           ! cm/yr/v_c,
                                                           ! v_c =
                                                           ! kappa/1000
                                                           ! (for km->m
                                                           ! conversion)
      end module sub_bench 
      
c *********************************************************************

      program conman
      use bc
      use fe
      use materials
      use solver
      use mstres
      use variables
c
      implicit double precision (a-h,o-z)
c
c.... remove above card for single-precision operation
c
        character*80 name1, name2, name3, name4, name5,
     &               name6, name7, name8, name9, name10
        
        include 'common.h'
c
c.... set up the input/output files
c
        read(5,1001) name1
        read(5,1001) name2
        read(5,1001) name3
        read(5,1001) name4
        read(5,1001) name5
        read(5,1001) name6
        read(5,1001) name7
        read(5,1001) name8
        read(5,1001) name9
        read(5,1001) name10

 1001   format( a80 )
        open (unit=iin,    file=name1 ,      status='old')
        open (unit=igeom,  file=name2 ,      status='old')
        open (unit=iout,   file=name3 ,      status='unknown')
        open (unit=irsin,  file=name4 ,      status='unknown')
        open (unit=irsout, file=name5 ,      status='unknown')
        open (unit=itsout, file=name6 ,      status='unknown') 
        open (unit=itout,  file=name7 ,      status='unknown')
        open (unit=imout,  file=name8 ,      status='unknown')
        open (unit=icomp,  file=name9 ,      status='unknown') 
        open (unit=igeoid, file=name10,      status='unknown')
c
c.... Input phase
c
      nmplot = 1    ! counter for movie files e.g. mov0001.ppm refer to 'ppmprt'
      call timer ('begin   ')
      call input
      close (iin)
      close (igeom)
c
c.... Execution phase
c
c.... check for execution mode
c
      if (iflow .ne. 0) then
c
c.... call the iterative solver--time integration driver routine
c
      call time_driver ( )
c
      call timer ('end     ')
c
      write(iout,2000) 
      write(iout,2001) ntotal
      rndpbyte=ntotal*8.e-6
      write(iout,2002) rndpbyte
c
      endif
c
 2000 format(
     & ' M e m o r y___ A l l o c a t i o n___ S t a t i s t i c s',//)
 2001 format('ntotal__= ',i12  ,'___Total_Words_of_Memory_Program ')
 2002 format('ndpbyte_= ',f12.2,'___Mb_Memory_Allocated')
c
c.... Close files
c
      close (iout)
      close (irsin)
      close (irsout)
      close (itsout)
      close (itout)
      close (imout)
      close (icomp)
      close (igeoid)
c
c.... end of the program
c  
      call exit (istat)
c
c.... set free of memories
c     
      deallocate (idv, vbc, idt, tbc, isuf, suf)
      deallocate (ien, lmv, lmt, shl, shdx, shdy, det)
      deallocate (dens, tref, mat, visc, diff, ra, alam, dmhu, tcon)
      deallocate (vlhsa, vlhsb, vrhs, vbcr, idiag, trhs)
      deallocate (tlhsa, tlhsb, idiagt)          ! Implicit and Picard
      deallocate (tlhs)                          ! Explicit
      deallocate (nb, stress, pmass)
      deallocate (x, v, t, tdot)
      end

c----------------------------------------------------------------------
      block data

      implicit double precision (a-h,o-z)

      common /io    / iin,    igeom, iout , itsout , itout , imout, 
     &                irsin , irsout, icomp, igeoid
c
c.... useful constants
      common /const / zero ,   pt25,     pt33,   pt5,    pt66,   one,
     &               onept5,   two ,    three,  four,  sixten,  eps7
c
      data       iin , igeom  , iout , itsout , itout , imout ,
     &         irsin , irsout , icomp, igeoid 
     
     &      /      10,     11 ,     12,     13,     14,     15,
     &             16,     17 ,     18,     19 /

c
      data        zero,                  pt25,                 
     &            pt33,                  pt5,                  
     &            pt66,                  one,                  
     &            onept5,                two,                  
     &            three,                 four,                 
     &            sixten,                eps7
     &      /     0.000000000000000D+0,   2.500000000000000D-1, 
     &            3.333333333333333D-1,   5.000000000000000D-1,   
     &            6.666666666666666D-1,   1.000000000000000D+0,   
     &            1.500000000000000D+0,   2.000000000000000D+0,   
     &            3.000000000000000D+0,   4.000000000000000D+0,   
     &            1.600000000000000D+1,   1.000000000000000D-7    /
      end
c
c----------------------------------------------------------------------
   
