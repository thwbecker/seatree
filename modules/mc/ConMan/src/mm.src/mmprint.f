*dk,mmprint
      subroutine mmprint()
C
C
C ######################################################################
C
C     PURPOSE:
C
C        This routine prints a map of the Memory Manger.
C
C     INPUT ARGUMENTS:
C
C        NONE
C
C     OUTPUT ARGUMENTS:
C
C        NONE
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
C        $Log:   /pvcs.config/utilities/src/mmprint.f_a  $
CPVCS    
CPVCS       Rev 1.7   Wed Jan 03 12:17:48 1996   het
CPVCS    Change the address calculation for UNICOS
CPVCS    
CPVCS       Rev 1.6   11/16/95 09:55:46   dcg
CPVCS    print type
CPVCS
CPVCS       Rev 1.5   11/09/95 16:55:06   dcg
CPVCS    add integer/real capability to memory manager
CPVCS
CPVCS       Rev 1.4   05/30/95 13:46:40   ejl
CPVCS    Cleaned up, Implicit none.
CPVCS
CPVCS
CPVCS       Rev 1.3   02/08/95 17:57:10   het
CPVCS    Change the table print format
CPVCS
CPVCS       Rev 1.2   01/26/95 11:26:04   het
CPVCS    Correct some character comparison problems.
CPVCS
CPVCS
CPVCS       Rev 1.1   01/04/95 21:59:16   llt
CPVCS    unicos changes (made by het)
CPVCS
CPVCS       Rev 1.0   11/10/94 12:43:28   pvcs
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
      integer i, j, k, numactive, icount, next, ierror
      integer ilast, iflag, iblkmin, jblkmin, kblkmin, iblkminsave
      character*32 blkin,prtin
C
      integer length, ipointer, itype
C
C ######################################################################
C
C
C
C     ******************************************************************
C     Get the maximum block address.
C
      iblkmin = 0
      numactive = 0
      do i=1,nblks
         if (mmsblk(5,i) .ge. 1) then
            iblkmin = max(iblkmin,mmsblk(3,i))
            numactive = numactive + 1
         endif
      enddo
      iblkminsave = iblkmin
C
C     ******************************************************************
C     Print the Memory Manager Data in index order.
C
      write(*,9000) numactive
 9000 format(/, 'Total number of active blocks: ',i10, /)
C
      write (*,9020)
 9020 format (/, 'INDEX', 1x, 'LENGTH', 1x, 'TYPE', 1x,'POINTER', 1x,
     &        'BLOCK_NAME', 21x, 'PARTITION_NAME', /)
C
      icount = 0
      do i=1,nblks
         if (mmsblk(5,i) .ge. 1) then
            icount = icount+mmsblk(4,i)*mmsblk(5,i)
            iblkmin = min(iblkmin,mmsblk(3,i))
            write(*,9040) i, mmsblk(4,i), mmsblk(5,i),mmsblk(3,i),
     &                       mmsnam(1,i), mmsnam(2,i)
 9040       format(1x, i4, i7, i2, 1x, i9, 2x, a32, a32)
         endif
      enddo
C
      write(*,9060) icount
 9060 format(/, 'Number of integer words used: ',i10, //)
C
C     ******************************************************************
C     Print the Memory Manager Data in address order.
C
      write(*,9000) numactive
C
      write (*,9020)
C
      ilast  = 0
      icount = 0
      ierror = 0
C
      do while ((icount.lt.numactive) .and. (ierror.eq.0))
C
         iflag = 0
         j=0
         do while ((j.lt.nblks) .and. (iflag.eq.0) .and. (ierror.eq.0))
            j=j+1
            if ((iblkmin.eq.mmsblk(3,j)) .and. (mmsblk(5,j).ge.1)) then
               iflag = j
               blkin    = mmsnam(1,j)
               prtin    = mmsnam(2,j)
               ipointer = mmsblk(3,j)
               length   = mmsblk(4,j)
               itype    = mmsblk(5,j)
               next = ipointer+NPTRFAC_INT*(mmsdim+itype*length+1)
               k=0
               do while ((k.lt.nblks) .and. (iflag.gt.0))
                  k=k+1
                  if ((k.ne.j) .and. (mmsblk(5,k).ge.1) .and.
     &                (mmsblk(3,k).gt.iblkmin)) then
                     if (mmsblk(3,k) .lt. next) then
                        iflag=0
                        write(*,9080) mmsnam(1,j), mmsnam(2,j),
     &                         ipointer,length,itype,
     &                                mmsnam(1,k), mmsnam(2,k),
     &                            mmsblk(3,k),mmsblk(4,k),mmsblk(5,k)
 9080                   format('The memory for ',a,1x,a /
     &                        i9,i7,i2/
     &                         ' overlaps the memory for ',a,1x,a/
     &                         i9,i7,i2)
                        ierror=1
                     endif
                  endif
               enddo
            endif
         enddo
C
         if ((iflag.eq.0) .and. (ierror.eq.0)) then
C
C....       There is a gap in the Memory Manager Storage.
C           Must find the next Block after the gap.
C
            jblkmin = iblkminsave
            j=0
            do while ((j.lt.nblks) .and. (ierror.eq.0))
               j=j+1
               if (mmsblk(5,j) .ge. 1) then
                  kblkmin = mmsblk(3,j)
                  length  = mmsblk(4,j)
                  itype   = mmsblk(5,j)
                  if(j.ne.ilast.and.iblkmin.gt.kblkmin.and.
     &              iblkmin.lt.(kblkmin+
     &                             NPTRFAC_INT*(mmsdim+itype*length+1)))
     &                   then
                     write(*,9080) mmsnam(1,j), mmsnam(2,j),
     &                         ipointer,length,itype,
     &                             mmsnam(1,ilast), mmsnam(2,ilast),
     &                            mmsblk(3,k),mmsblk(4,k),mmsblk(5,k)
                     ierror=2
                  elseif(mmsblk(3,j).gt.iblkmin) then
                     iflag = j
                     blkin    = mmsnam(1,j)
                     prtin    = mmsnam(2,j)
                     ipointer = mmsblk(3,j)
                     length   = mmsblk(4,j)
                     itype    = mmsblk(5,j)
                     jblkmin  = min(jblkmin,mmsblk(3,j))
                  endif
               endif
            enddo
C
            iblkmin = jblkmin
C
         else
C
            if (iflag .gt. 0) then
               icount = icount+1
               write(*,9040) iflag, length, itype,ipointer, blkin, prtin
               ilast = iflag
               iblkmin = ipointer+NPTRFAC_INT*(mmsdim+itype*length+1)
            endif
C
         endif
C
      enddo
C
      if(ierror.gt.0) then
         write(*,'(a)') 'Error in memory map: Last good block is'
         write(*,9040) iflag, length, itype,ipointer, blkin, prtin
      endif
C
      return
      end
