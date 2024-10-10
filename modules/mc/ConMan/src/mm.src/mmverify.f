*dk,mmverify
      subroutine mmverify()
C
C
C ######################################################################
C
C     PURPOSE -
C
C        Verifies the Memory Management Storage.
C
C     INPUT ARGUMENTS -
C
C        NONE
C
C     OUTPUT ARGUMENTS -
C
C        NONE
C
C     CHANGE HISTORY -
C
C        $Log:   /pvcs.config/utilities/src/mmverify.f_a  $
CPVCS    
CPVCS       Rev 1.4   Wed Jan 03 12:17:54 1996   het
CPVCS    Change the address calculation for UNICOS
CPVCS    
CPVCS       Rev 1.3   11/21/95 13:24:44   dcg
CPVCS    allow for real and integer types -- mmsblk is integer
CPVCS
CPVCS       Rev 1.2   05/30/95 13:47:08   ejl
CPVCS    Cleaned up, Implicit none.
CPVCS
CPVCS
CPVCS       Rev 1.1   01/04/95 21:59:32   llt
CPVCS    unicos changes (made by het)
CPVCS
CPVCS       Rev 1.0   11/10/94 12:43:42   pvcs
CPVCS    Original version.
C
C ######################################################################
C
      implicit none
C
C ######################################################################
C
 
C$Log:   /pvcs.config/utilities/src/machine.h_a  $
CPVCS    
CPVCS       Rev 1.6   11/28/95 11:39:12   dcg
CPVCS    add NBYTES_real NBYTES_int for cray/workstation
CPVCS    implementation
CPVCS
CPVCS       Rev 1.3   11/27/95 11:36:38   het
CPVCS    Add the NPTRFAC_INT parameter to indicate the length of integers.
CPVCS
CPVCS       Rev 1.4   11/27/95 11:32:56   het
CPVCS    Add the NPTRFAC_INT parameter to indicate the length of integers
CPVCS
CPVCS       Rev 1.3   06/06/95 16:08:32   dcg
CPVCS
CPVCS       Rev 1.2   06/06/95 16:01:54   dcg
CPVCS    add type statments
C                           comdeck machine
 
 
C                 PARAMETERize Cray-Machine  and CTSS
C                     System-Dependent Quantities.
 
 
C                       CRAY Machine Quantities
 
C     KNCPW    The number of characters per computer (REAL*8) word.
      INTEGER KNCPW
               PARAMETER (KNCPW = 8)
 
C     KNWPN    The number of computer (REAL*8) words per character name.
      INTEGER KNWPN
               PARAMETER (KNWPN = 4)
 
C     KNCPN    The number of characters per character name.
      INTEGER KNCPN
               PARAMETER (KNCPN = KNCPW * KNWPN)
 
C     KNBPW    The number of bits per computer (REAL*8) word.
      INTEGER KNBPW
               PARAMETER (KNBPW = 64)
 
 
C                       CTSS System  Quantities
 
C     KWLSECT  The length in computer words of a disk sector.
      INTEGER KWLSECT
               PARAMETER (KWLSECT = 512)
 
C
C     nptrfac The multiplier for pointer arithmetic.  i.e for:
C             WORKSTATION set NPTRFAC=8, NPTRFAC_INT=4
C             CRAY        set NPTRFAC=1, NPTRFAC_INT=1
      INTEGER NPTRFAC, NPTRFAC_INT
              parameter (NPTRFAC=8, NPTRFAC_INT=4)
C
C     NBYTES_REAL number of bytes in a real variable
C     NBYTES_INT number of bytes in a integer variable
C             WORKSTATION set NBYTES_REAL = 8, NBYTES_INT = 4
C             CRAY set NBYTES_REAL = 8, NBYTES_INT = 8
      integer NBYTES_REAL, NBYTES_INT
              parameter (NBYTES_REAL = 8, NBYTES_INT=4)
C                         ( end  of  machine )
*cd,mmsl
C
C
C#######################################################################
C
C     PURPOSE:
C
C        Common Deck for Memory Manager Routines.
C
C     CHANGE HISTORY -
C
C        $Log:   /pvcs.config/utilities/src/mmsl.h_a  $
CPVCS    
CPVCS       Rev 1.7   11/29/95 10:01:58   dcg
CPVCS    use NBYTES_REAL and NBYTES_INT for nwlen2, nwlen3
CPVCS
CPVCS       Rev 1.6   11/27/95 13:17:32   het
CPVCS    Add the NPTRFAC_INT parameter for integer length assignments.
CPVCS
CPVCS       Rev 1.5   11/14/95 16:59:18   dcg
CPVCS    make header length 8
CPVCS
CPVCS       Rev 1.4   05/30/95 13:57:40   ejl
CPVCS    Cleaned up, Implicit none.
CPVCS
CPVCS
CPVCS       Rev 1.3   05/24/95 15:07:12   ejl
CPVCS    Typed all variables.
CPVCS
C
C#######################################################################
C
c     include "machine.h"
C
C#######################################################################
C
C
C     FOR A CRAY SET nwlen2 = 1
C     FOR A UNIX WORKSTATION SET nwlen2=8
C
C    SET nwlen1 to length in bytes of real data
C    SET nwlen3 to length in byts of integer data
C
C   mmsdim must always be even  because of packing problems and
C   work alignment considerations
C
      integer mmsdim, nwlen1, nwlen2, nwchar1, nwlen3
      parameter(mmsdim=8, nwlen1=8,
     *          nwlen2=NBYTES_REAL, nwlen3=NBYTES_INT,
     *          nwchar1=32)
C
      common /mmsdat/maxblks, nblks, mmsindptr, mmsnamptr
      integer        maxblks, nblks
C
      pointer(mmsindptr, mmsblk)
      integer mmsblk(mmsdim,*)
C
      pointer(mmsnamptr, mmsnam)
      character*32 mmsnam(2,*)
C
C     ------------------------------------------------------------------
C     The array mmsblk/mmsnam is the memory manager control block and
C     stores the following data:
C         mmsnam(1,i) - block name.
C         mmsnam(2,i) - partition name.
C         mmsblk(3,i) - block address.
C         mmsblk(4,i) - number of words allocated.
C         mmsblk(5,i) - delete/type flag (0=deleted, 1=integer, 2=real*8).
C         mmsblk(6,i) - unused included for word alignment
C         mmsblk(7,i) - unused included for word alignment
C         mmsblk(8,i) - unused included for word alignment
C     ------------------------------------------------------------------
C
 
C
C ######################################################################
C
      integer i, j, iflagg, iptruser, length, istatus, iflag1, iflag2
      integer itype
C
      character*32 blkout, prtout
C
      pointer (ipointer, xarray)
      real*8 xarray(mmsdim+1)
      pointer (ipointer, iarray)
      integer iarray(mmsdim+1)
C
      integer idebug
      data idebug / 1 /
C
C ######################################################################
C
C
C
C     ******************************************************************
C     Check that this block exists.
C
      if(idebug.ne.0) then
C
         iflagg=0
C
         do i=1,nblks
C
            blkout   = mmsnam(1,i)
            prtout   = mmsnam(2,i)
            iptruser = mmsblk(3,i)
            length   = mmsblk(4,i)
            itype    = mmsblk(5,i)
            istatus  = mmsblk(5,i)
C
            ipointer=iptruser-NPTRFAC_INT*mmsdim
C
            iflag1=0
            iflag2=0
C
            if(istatus.ge.1) then
C
C              *********************************************************
C              Active Block, check header.
C
               do j=1,mmsdim
                  if(mmsblk(j,i).ne.iarray(mmsdim-j+1)) then
                     iflag1=iflag1+1
                     print *,'Memory management error: Block=',blkout,
     *                       ' Part=',prtout,' Length=',length,
     *                       ' Address=',iptruser
                     print *,'   Header was overwritten'
                  endif
               enddo
C
C              *********************************************************
C              Active Block, check trailer.
C
               if(mmsblk(1,i).ne.iarray(mmsdim+length*itype+1)) then
                  iflag2=1
                  print *,'Memory management error: Block=',blkout,
     *                    ' Part=',prtout,' Length=',length,
     *                    ' Address=',iptruser
                  print *,'   Trailer was overwritten'
                  print *,' Block Number = ',i,' mmsblk ',mmsblk(1,i)
                  print *,' iarray ',iarray(mmsdim+itype*length+1)
               endif
C
            endif
C
            iflagg=max(iflagg,iflag1,iflag2)
C
         enddo
C
         if(iflagg.ne.0) then
            call mmprint()
            call heapexit()
         endif
C
      endif
C
      return
      end
