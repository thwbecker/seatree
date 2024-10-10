*dk,mmrelblk
      subroutine mmrelblk(blkin,prtin,iadr,icscode)
C
C
C ######################################################################
C
C     PURPOSE:
C
C        This routine releases memory previously allocated for a block.
C        After the memory is released, the block is deleted from the
C        memory manager control block sets.
C
C     INPUT ARGUMENTS:
C
C        blkin - the block's name.
C        prtin - the block's partition name.
C        len   - the block's new length in computer words.
C        iadr  - a possible address of the block, not used.
C
C     OUTPUT ARGUMENTS:
C
C        icscode - completion status code:
C                    0 - normal
C                  -16 - block does not exists.
C                   -5 - system memory manager error.
C                   -6 - memory manager control block error.
C
C     COMMON VARIABLES:
C
C        maxblks - no. of control block sets allocated.
C        nblks - no. of control block sets used.
C        mmsindptr - pointer to the control block sets.
C        mmsnamptr - pointer to the control block sets.
C        mmsdim - no. of control words in a control block set.
C
C     CHANGE HISTORY -
C
C        $Log:   /pvcs.config/utilities/src/mmrelblk.f_a  $
CPVCS    
CPVCS       Rev 1.7   Wed Jan 03 12:17:52 1996   het
CPVCS    Change the address calculation for UNICOS
CPVCS    
CPVCS       Rev 1.6   11/09/95 16:55:08   dcg
CPVCS    add integer/real capability to memory manager
CPVCS
CPVCS       Rev 1.5   05/30/95 13:46:42   ejl
CPVCS    Cleaned up, Implicit none.
CPVCS
CPVCS
CPVCS       Rev 1.4   05/24/95 14:56:52   ejl
CPVCS    Changed the way charter length was done.
CPVCS
CPVCS
CPVCS       Rev 1.3   05/15/95 16:08:36   dcg
CPVCS    fix character length problem
CPVCS
CPVCS       Rev 1.2   03/17/95 21:16:18   het
CPVCS    Changes supporting the RDBMS/CAS systems
CPVCS
CPVCS       Rev 1.1   01/04/95 21:59:20   llt
CPVCS    unicos changes (made by het)
CPVCS
CPVCS       Rev 1.0   11/10/94 12:43:30   pvcs
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
      character*(*) blkin, prtin
      integer iadr
C
      integer icscode
C
C ######################################################################
C
      integer i, ifound, length, itype
      integer ierror
C
      pointer (ipointer, xarray)
      real*8 xarray(mmsdim+1)
      pointer (ipointer, iarray)
      integer iarray(mmsdim+1)
C
C ######################################################################
C
      integer freef
C
C ######################################################################
C
C
C
C     ******************************************************************
C     Verify the integrity of the memory before doing anything.
      if(nblks.gt.0) call mmverify()
C
C     ******************************************************************
C     Get the Index to the Block.
C
      call mmgetindex(blkin,prtin,ifound,icscode)
C
      if (ifound .eq. 0) then
C
C        ***************************************************************
C        Block does not exist.
C
         icscode=-16
         print 9001, blkin, prtin
 9001    format('  Error, mmrelblk - block ',a,1x,a,
     &          ' does not exist')
C
      elseif (mmsblk(5,ifound).eq.0) then
C
C        ***************************************************************
C        Inactive block:
C
         icscode=-21
         print 9002, blkin, prtin
 9002    format('  Error, mmrelblk - trying to release an inactive',
     *          ' block',a,1x,a)
         call heapexit()
C
      else
C
C        ***************************************************************
C        Get the pointer, free memory and set delete flag.
C
         ipointer = mmsblk(3,ifound)-NPTRFAC_INT*mmsdim
         length   = mmsblk(4,ifound)
         itype    = mmsblk(5,ifound)
C
         do i=1,mmsdim+length+1
            iarray(i)=0
         enddo
C
C....    Free the memory.
C
         ierror=freef(ipointer)
C
C....    Make the Block Inactive.
C
         mmsblk(5,ifound)=0
C
C        ***************************************************************
C        Compress the memory manager control block.
C
         call mmcomprs(0,icscode)
C
      endif
C
      return
      end
