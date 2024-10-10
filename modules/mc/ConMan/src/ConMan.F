      program conman

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
c **             Incompressible Thermal Convective Flows.            **
c **                                                                 **
c **                                                                 **
c ** written by:     Arthur Raefsky and Scott King .                 **
c **                                                                 **
c **                                                                 **
c *********************************************************************
c *********************************************************************
c----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
c.... remove above card for single-precision operation
c
        character*80 name1,name2,name3,name4,name6,name7,name8,name9
        common /temp1 / x(50000)
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
 1001   format( a80 )
        open (unit=iin,    file=name1 , status='old')
        open (unit=igeom,  file=name2 , status='old')
        open (unit=iout,   file=name3 , status='unknown')
        open (unit=irsin,  file=name4 , status='unknown')
c       open (unit=irsout, file=name5 , status='unknown')
        open (unit=itsout, file=name6 , status='unknown' ) 
        open (unit=itout,  file=name7 , status='unknown')
        open (unit=imout,  file=name8 , status='unknown' ) 
        open (unit=igeoid,  file=name9 , status='unknown' ) 
c
c.... Input phase
c
      call timer ('begin   ')
      call input
      close (iin)
      close (igeom)
c
c.... Execution phase
c
      call proces
      call timer ('end     ')

      write(iout,2000) 
      write(iout,2001) ntotal
      rndpbyte=ntotal*4.e-6
      write(iout,2002) rndpbyte

      call mmprint(iout)

 2000 format(
     & ' M e m o r y   A l l o c a t i o n   S t a t i s t i c s',//)
 2001 format('ntotal   = ',i12  ,'   Total Words of Memory Program ')
 2002 format('ndpbyte  = ',f12.2,'   Mb Memory Allocated')

c
c.... Close files
c
      close (iout)
      close (irsin)
      close (itsout)
      close (itout)
      close (imout)
c
c.... end of the program
c
      call exit (istat)
      end
      block data
c----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
c.... remove above card for single-precision operation
c
c.... dynamic storage pointer management data
c
c     common /point / mbeg,   mend,   mprec,  mdct,   ndict
c
c.... io channels
      common /io    / iin  ,    igeom,  iout  , itsout , itout ,
     &                imout ,   irsin , irsout, igeoid
c
c.... useful constants
      common /const / zero ,   pt25,     pt33,   pt5,    pt66,   one,
     &               onept5,   two ,    three,  four,  sixten,  eps7
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c
c     data      mbeg,    mend,    mprec,   mdct,    ndict
c    &      /      1,130000000,        2,      1,        0 /
c
      data       iin , igeom  , iout , itsout , itout , imout ,
     &         irsin , irsout , igeoid
     &      /      10,     11 ,     12,     13,     14,     15,
     &             16,     17 ,     18/
c
      data        zero,                  pt25,                 
     &            pt33,                  pt5,                  
     &            pt66,                  one,                  
     &            onept5,                two,                  
     &            three,                 four,                 
     &            sixten,                eps7
     &      /     0.00000000000000d+0,   2.50000000000000d-1, 
     &            3.33333333333333d-1,   5.00000000000000d-1,   
     &            6.66666666666666d-1,   1.00000000000000d+0,   
     &            1.50000000000000d+0,   2.00000000000000d+0,   
     &            3.00000000000000d+0,   4.00000000000000d+0,   
     &            1.60000000000000d+1,   1.00000000000000d-7    /
c
c----------------------------------------------------------------------
      end