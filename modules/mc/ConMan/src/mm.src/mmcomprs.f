*dk,mmcomprs
      subroutine mmcomprs(iholesize,icscode)
C
C
C ######################################################################
C
C     PURPOSE:
C
C        This routine removes memory manager control blocks marked as
C        deleted and compress the memory manager control block sets
C
C     INPUT ARGUMENTS:
C
C        none
C
C     OUTPUT ARGUMENTS:
C
C        icscode - completion status code:
C                    0 - normal.
C
C     COMMON VARIABLES:
C
C        maxblks   - no. of control block sets allocated.
C        nblks     - no. of control block sets used.
C        mmsindptr - pointer to the control block sets.
C        mmsnamptr - pointer to the control block sets.
C        mmsdim    - no. of control words in a control block set.
C
C     CHANGE HISTORY -
C
C        $Log:   /pvcs.config/utilities/src/mmcomprs.f_a  $
CPVCS    
CPVCS       Rev 1.6   11/15/95 09:36:44   dcg
CPVCS    use loop on mmsdim instead of individual indices
CPVCS
CPVCS       Rev 1.5   11/14/95 16:59:26   dcg
CPVCS    make header length 8
CPVCS
CPVCS       Rev 1.3   06/28/95 11:10:16   het
CPVCS    Correct an error when zero blocks are compressed
CPVCS
CPVCS       Rev 1.2   05/30/95 13:46:08   ejl
CPVCS    Cleaned up, Implicit none.
CPVCS
CPVCS
CPVCS       Rev 1.1   01/04/95 21:58:44   llt
CPVCS    unicos changes (made by het)
CPVCS
CPVCS       Rev 1.0   11/10/94 12:43:04   pvcs
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
C# #####################################################################
C
      integer iholesize
      integer icscode
C
C ######################################################################
C
      integer i, j, ndel, istrt, iend, irepl, k
C
C ######################################################################
C
C
C
      icscode=0
C
C     ******************************************************************
C     Verify the integrity of the memory before doing anything.
      if(nblks.gt.0) call mmverify()
C
C     ******************************************************************
C     Find first deleted block and count deleted blocks.
C
      ndel=0
      istrt=0
      do i=1,nblks
         if (mmsblk(5,i) .eq. 0) then
            ndel=ndel+1
            if (istrt.eq.0) istrt=i
         endif
      enddo
C
C     ******************************************************************
C     If there are not blocks to delete then bailout.
C
      if(ndel.gt.0) then
C
C
C        ***************************************************************
C        Replace deleted blocks with non-deleted blocks.
C
         iend=nblks-ndel
C
         do i=istrt,iend
C
C           ------------------------------------------------------------
C           See if this is a deleted block.
C
            if (mmsblk(5,i) .eq. 0) then
C
C              .........................................................
C              Search backwards to find next non-deleted block.
C
               irepl=0
               j=nblks+1
               do while ((j.gt.i) .and. (irepl.eq.0))
                  j=j-1
                  if (mmsblk(5,j) .ge. 1) irepl=j
               enddo
C
C              .........................................................
C              Replace deleted block with next non-deleted block.
C
               if (irepl .gt. 0) then
C
                  mmsnam(1,i) = mmsnam(1,irepl)
                  mmsnam(2,i) = mmsnam(2,irepl)
                  do k=1,mmsdim
                     mmsblk(k,i) = mmsblk(k,irepl)
                  enddo
C
C                 ......................................................
C                 Mark moved block as deleted.
C
                  mmsnam(1,irepl) = ' '
                  mmsnam(2,irepl) = ' '
                  mmsblk(1,irepl) = 0
                  mmsblk(2,irepl) = 0
                  mmsblk(3,irepl) = 0
                  mmsblk(4,irepl) = 0
                  mmsblk(5,irepl) = 0
C
                  nblks = irepl-1
C
               endif
C
            endif
C
         enddo
C
C        ...............................................................
C        Reset the number of bloks.
C
         nblks=iend
C
      endif
C
C
      return
      end
