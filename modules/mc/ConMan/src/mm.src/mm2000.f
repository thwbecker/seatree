!!
!!  Memory Manager 2000 (:-P)
!!
!!  Neil Carlson <Neil.Carlson@Motorola.com>
!!  Motorola Computational Materials Group @ LANL
!!  6 April 2000
!!
!!  This package is intended to be a drop-in replacement for the present
!!  memory management routines (mm* files) used by LaGriT.  It maintains
!!  the interface and aims to replicate the behavior of the previous
!!  routines, while reimplementing the internals and clearing away years
!!  of accumulated cruft.
!!
!!
!!  PROCEDURES
!!
!!  The following procedure arguments are used throughout.
!!
!!    character*(*) aname -- Name of the managed array
!!    character*(*) pname -- Partition name of the managed array
!!    pointer (aptr, array(*)) -- Cray-style pointer to managed array
!!    integer length -- Size of the managed array
!!    integer type -- Type of the managed array:
!!                    1 == integer
!!                    2 == real*8
!!                    3 == character*32
!!    integer status -- Return status variable.  STATUS = 0 indicates
!!       a successful operation.
!!
!!  1. PUBLIC ALLOCATION PROCEDURES
!!
!!  MMGETBLK (ANAME, PNAME, APTR, LENGTH, TYPE, STATUS)
!!    Allocates a new managed array of type TYPE and size LENGTH.
!!    Input: ANAME, PNAME, LENGTH, TYPE
!!    Output: APTR, STATUS
!!
!!  MMNEWLEN (ANAME, PNAME, APTR, LENGTH, STATUS)
!!    Reallocates an existing managed array to size LENGTH.
!!    Input: ANAME, PNAME, LENGTH
!!    Output: APTR, STATUS
!!
!!  MMRELBLK (ANAME, PNAME, APTR, STATUS)
!!    Deallocates a managed array.
!!    Input: ANAME, PNAME
!!    Output: STATUS
!!    APTR is unused.
!!
!!  MMRELPRT (PNAME, STATUS)
!!    Deallocates all managed arrays in a partition.
!!    Input: PNAME; Output: STATUS
!!
!!  MMNAMPRT (NEWPNAME, OLDPNAME, STATUS)
!!    Changes the name of partition OLDPNAME to NEWPNAME.
!!    Input: NEWPNAME, OLDPNAME.  Output: STATUS.
!!
!!  MMINCBLK (ANAME, PNAME, APTR, INCR, STATUS)
!!    Reallocates an existing managed array incrementing its size
!!    by INCR.
!!    Input: ANAME, PNAME, INCR
!!    Output: APTR, STATUS
!!
!!  MMGGETBK (ANAME, PNAME, APTR, LENGTH, TYPE, STATUS)
!!    If the named array does not exist, allocate a new managed array
!!    of type TYPE and size LENGTH.  If the named array exists and
!!    LENGTH exceeds its current size, reallocate it with size LENGTH.
!!    Otherwise simply return the current array pointer.
!!    Input: ANAME, PNAME, LENGTH, TYPE
!!    Output: APTR, STATUS
!!
!!  MM_OVALL (ANAME, PNAME, APTR, NEED, EXTRA, LENGTH, TYPE, STATUS)
!!    If LENGTH < NEED and LENGTH = 0, a new managed array is allocated
!!    of type TYPE and size NEED + EXTRA.  IF LENGTH < NEED and
!!    LENGTH > 0, the existing named array is reallocated with size
!!    NEED + EXTRA.  In either case the new size is returned in LENGTH.
!!    For proper operation, APTR must point to the named array (if it
!!    exists) and LENGTH must be its size; LENGTH = 0 indicates the
!!    named array does not exist.
!!    Input: ANAME, PNAME, NEED, EXTRA, LENGTH, TYPE
!!    Output: APTR, LENGTH, STATUS
!!
!!  MMCMPRS (iholesize, status)
!!    This does absolutely nothing but return STATUS = 0.  The previous
!!    version compressed the array descriptor blocks, but in this new
!!    implementation this operation is moot.
!!
!!  2. PUBLIC INQUIRY PROCEDURES
!!
!!  MMBLKLEN (ANAME, PNAME, APTR, LENGTH, STATUS)
!!  MMFINDBK (ANAME, PNAME, APTR, LENGTH, STATUS)
!!    Returns a pointer to, and size of, the named managed array.
!!    Input: ANAME, PNAME.  Output: APTR, LENGTH
!!
!!  MMGETINDEX (ANAME, PNAME, INDEX, STATUS)
!!    Returns the index of the named managed array.  INDEX > 0 indicates
!!    the array exists, INDEX = 0 it does not.
!!    Input: ANAME, PNAME.  Output: INDEX, STATUS
!!
!!  MMGETLEN (APTR, LENGTH, STATUS)
!!    Returns the size of a managed array by address.
!!    Input: APTR.  Output: LENGTH, STATUS.
!!
!!  MMGETNAM (APTR, ANAME, PNAME, STATUS)
!!    Returns the name and partition of a managed array by address.
!!    Input: APTR.  Output: ANAME, PNAME, STATUS
!!
!!  MMGETPR (ANAME, PNAME, APTR, STATUS)
!!    Returns the pointer to a named managed array.
!!    Input: ANAME, PNAME.  Output: APTR, STATUS
!!
!!  MMGETTYP (APTR, TYPE, STATUS)
!!    Returns the type of a managed array by address.
!!    Input: APTR.  Output: TYPE, STATUS
!!
!!  3. DIAGNOSTIC PROCEDURES
!!
!!  MMVERIFY ()
!!    This aims to verify the integrity of the managed array storage
!!    by checking that the prolog and epilog blocks have not been
!!    overwritten.  The code will print out an array map and halt if
!!    it detects such memory corruption.
!!
!!  MMPRINT ()
!!    Print out an address map of the managed arrays.
!!
!!  4. PRIVATE PROCEDURES -- DO NOT USE
!!
!!  MM_NEW_AD (INDEX, STATUS)
!!    Return an index (in the AD array) for a new array descriptor.
!!
!!  MM_INDEX_BY_NAME (ANAME, PNAME, INDEX)
!!    Return the AD-array index of the named managed array.
!!
!!  MM_INDEX_BY_ADDR (APTR, INDEX)
!!    Return the AD-array index of the managed array by address.
!!
!!  MM_PANIC (PROC, MESG)
!!    Print error message and halt.
!!
!!  MM_WARN (PROC, MESG)
!!    Print warning message and return.
!!
!!
!!  CHANGE HISTORY
!!
!!  $Log:   /pvcs.config/utilities/src/mm2000.f_a  $
CPVCS    
CPVCS       Rev 1.2   Fri Apr 07 08:34:44 2000   dcg
CPVCS    rename blockdata
CPVCS    replace f90 syntax with f77 (specifically f77 does not recognize plain 'do')
CPVCS
CPVCS       Rev 1.1   Thu Apr 06 13:15:40 2000   nnc
CPVCS    Fixed typo in mm_ovall name.
CPVCS
CPVCS       Rev 1.0   Thu Apr 06 10:49:48 2000   nnc
CPVCS    Initial revision.
 
 
      block data mm2000data
      integer max_ad, num_ad, first, next, adnptr, adaptr, lptr
      common /mm2000/ max_ad, num_ad, first, next, adnptr, adaptr, lptr
      data max_ad, num_ad, first, next, adnptr, adaptr, lptr /7*0/
      end
 
!!
!!  1. PUBLIC ALLOCATION PROCEDURES
!!
 
      subroutine mmgetblk (aname, pname, aptr, length, type, status)
 
      implicit none
 
      character*(*) aname, pname
      integer length, type, status
 
      pointer (aptr,int_array), (aptr,real_array), (aptr,char_array)
      integer int_array(*)
      real*8 real_array(*)
      character*32 char_array(*)
 
      pointer (pptr,prolog), (eptr,epilog)
      integer prolog(*), epilog(*)
 
      include 'mm2000.h'
 
      integer array_bytes, block_bytes, index, j
 
      integer mallocf
      external mallocf
 
      call mmverify ()
 
      call mm_index_by_name (aname, pname, index)
      if (index .gt. 0) then  ! Array already exists
        status = -17
        call mm_panic ('MMGETBLK',
     &    'Array ' // aname // ' in ' // pname // 'already exists')
      end if
 
      if (length .lt. 1) then
        status = -21
        call mm_panic ('MMGETBLK', 'Invalid length request')
      end if
 
      if (type .lt. 1 .or. type .gt. 3) then
        status = -21
        call mm_panic ('MMGETBLK', 'Invalid type request')
      end if
 
      !! Allocate memory for the array
      if (type .eq. 1) then
        array_bytes = BYTES_PER_INT * length
      elseif (type .eq. 2) then
        array_bytes = BYTES_PER_REAL * length
      elseif (type .eq. 3) then
        array_bytes = BYTES_PER_CHAR * length
      endif
      block_bytes = PROLOG_BYTES + array_bytes + EPILOG_BYTES
 
      pptr = mallocf(block_bytes)
      aptr = pptr + PROLOG_BYTES
      eptr = aptr + array_bytes
 
c changed sdk
      if (pptr .eq. 0) then ! Failed to allocate memory
        status = -5
        call mm_panic ('MMGETBLK', 'Unable to allocate memory')
      end if
 
      !! Get new array descriptor and initialize
      call mm_new_ad (index, status)
      ad_name(1,index) = aname
      ad_name(2,index) = pname
      ad_addr(1,index) = pptr
      ad_addr(2,index) = aptr
      ad_addr(3,index) = eptr
      ad_addr(4,index) = type
      ad_addr(5,index) = length
      ad_addr(6,index) = block_bytes
 
      !! Initialize the array
      if (type .eq. 1) then
        do j = 1, length
          int_array(j) = 0
        end do
      else if (type .eq. 2) then
        do j = 1, length
          real_array(j) = 0.0
        end do
      else if (type .eq. 3) then
        do j = 1, length
          char_array(j) = ' '
        end do
      end if
 
      !! Install the prolog and epilog; arbitrary, but specified, data.
      do j = 1, AD_SIZE
        prolog(j) = ad_addr(j,index)
        epilog(j) = ad_addr(j,index)
      end do
 
      status = 0
 
      return
      end
 
 
      subroutine mmnewlen (aname, pname, aptr, length, status)
 
      implicit none
 
      character*(*) aname, pname
      integer length, status
 
      pointer (aptr,int_array), (aptr,real_array), (aptr,char_array)
      integer int_array(*)
      real*8 real_array(*)
      character*32 char_array(*)
 
      pointer (pptr,prolog), (eptr,epilog)
      integer prolog(*), epilog(*)
 
      include 'mm2000.h'
 
      integer j, orig_length, type, array_bytes, block_bytes, index
 
      integer reallocf
      external reallocf
 
      status = 0
 
      call mmverify ()
 
      call mm_index_by_name (aname, pname, index)
      if (index .eq. 0) then  ! Array doesn't exist
        status = -16
        call mm_panic ('MMNEWLEN',
     &    'Array ' // aname // ' in ' // pname // 'does not exist')
      end if
 
      if (length .lt. 1) then
        status = -21
        call mm_panic ('MMNEWLEN', 'Invalid length request')
      end if
 
      pptr = ad_addr(1,index)
      aptr = ad_addr(2,index)
      type = ad_addr(4,index)
      orig_length = ad_addr(5,index)
 
      if (length .eq. orig_length) return
 
      !! Allocate memory for the array
      if (type .eq. 1) then
        array_bytes = BYTES_PER_INT * length
      elseif (type .eq. 2) then
        array_bytes = BYTES_PER_REAL * length
      elseif (type .eq. 3) then
        array_bytes = BYTES_PER_CHAR * length
      endif
      block_bytes = PROLOG_BYTES + array_bytes + EPILOG_BYTES
 
      pptr = reallocf(pptr, block_bytes)
      aptr = pptr + PROLOG_BYTES
      eptr = aptr + array_bytes
 
c changed sdk
      if (pptr .eq. 0) then ! Failed to allocate memory
        status = -5
        call mm_panic ('MMNEWLEN', 'Unable to allocate memory')
      end if
 
      !! Reinitialize array descriptor
      ad_addr(1,index) = pptr
      ad_addr(2,index) = aptr
      ad_addr(3,index) = eptr
      ad_addr(4,index) = type
      ad_addr(5,index) = length
      ad_addr(6,index) = block_bytes
 
      !! Initialize new part of array
      if (length .gt. orig_length) then
        if (type .eq. 1) then
          do j = orig_length + 1, length
            int_array(j) = 0
          end do
        else if (type .eq. 2) then
          do j = orig_length + 1, length
            real_array(j) = 0.0
          end do
        else if (type .eq. 3) then
          do j = orig_length + 1, length
            char_array(j) = ' '
          end do
        end if
      end if
 
      !! Reinstall the prolog and epilog
      do j = 1, AD_SIZE
        prolog(j) = ad_addr(j,index)
        epilog(j) = ad_addr(j,index)
      end do
 
      status = 0
 
      return
      end
 
 
      subroutine mmrelblk (aname, pname, aptr, status)
 
      implicit none
 
      character*(*) aname, pname
      integer status
 
      pointer (aptr, array)
      integer array(*)
 
      include 'mm2000.h'
 
      pointer (pptr, block)
      integer block(*)
 
      integer prev, index, j, freef ,i
 
      external freef
 
      call mmverify ()
 
      prev = 0
      index = first
      do i=1,100000000
        if (index .lt. 1) then
          status = -16
          call mm_warn ('MMRELBLK',
     &      'Array ' // aname // ' in ' // pname // ' does not exist')
          return
        end if
        if (aname .eq. ad_name(1,index) .and.
     &      pname .eq. ad_name(2,index)) then
 
          !! DEALLOCATE THE MEMORY
          pptr = ad_addr(1,index)
          status = freef(pptr)
 
          !! CLEAR THE ARRAY DESCRIPTOR
          ad_name(1,index) = ' '
          ad_name(2,index) = ' '
          do j = 1, AD_SIZE
            ad_addr(j,index) = 0
          end do
 
          !! MOVE THE ARRAY DESCRIPTOR TO THE FREE LIST
          if (prev .gt. 0) then
            link(prev) = link(index)
            link(index) = next
            next = index
          else
            first = link(index)
            link(index) = next
            next = index
          end if
          num_ad = num_ad - 1
          status = 0
          return
        end if
        prev = index
        index = link(index)
      end do
 
      return
      end
 
 
      subroutine mmrelprt (pname, status)
 
      implicit none
 
      character*(*) pname
      integer status
 
      include 'mm2000.h'
 
      pointer (pptr, block)
      integer block(*)
 
      integer prev, index, j, freef ,i
 
      external freef
 
      call mmverify ()
 
      prev = 0
      index = first
      status = -16
      do i=1,100000000
        if (index .lt. 1) then
          if (status .ne. 0) then
!            call mm_warn ('MMRELPRT',
!     &                    'Partition '// pname // ' does not exist')
          end if
          return
        end if
        if (pname .eq. ad_name(2,index)) then
 
          !! DEALLOCATE THE MEMORY
          pptr = ad_addr(1,index)
          status = freef(pptr)
 
          !! CLEAR THE ARRAY DESCRIPTOR
          ad_name(1,index) = ' '
          ad_name(2,index) = ' '
          do j = 1, AD_SIZE
            ad_addr(j,index) = 0
          end do
 
          !! MOVE THE ARRAY DESCRIPTOR TO THE FREE LIST
          if (prev .gt. 0) then
            link(prev) = link(index)
            link(index) = next
            next = index
            index = link(prev)
          else
            first = link(index)
            link(index) = next
            next = index
            index = first
          end if
          num_ad = num_ad - 1
          status = 0
        else
          prev = index
          index = link(index)
        end if
      end do
 
      return
      end
 
 
      subroutine mmnamprt (newpname, oldpname, status)
 
      implicit none
 
      character*(*) newpname, oldpname
      integer status
 
      include 'mm2000.h'
 
      integer j ,i
 
      j = first
      status = -16
      do i=1,100000000
        if (j .lt. 1) return
        if (oldpname .eq. ad_name(2,j)) then
          ad_name(2,j) = newpname
          status = 0
        end if
        j = link(j)
      end do
 
      return
      end
 
 
      subroutine mmincblk (aname, pname, aptr, incr, status)
 
      implicit none
 
      character*(*) aname, pname
      integer incr, status
 
      pointer (aptr, array)
      integer array(*)
 
      integer length
 
      call mmfindbk (aname, pname, aptr, length, status)
      if (status .eq. 0) then
        length = length + incr
        call mmnewlen(aname, pname, aptr, length, status)
      end if
 
      return
      end
 
 
      subroutine mmggetbk (aname, pname, aptr, length, type, status)
 
      implicit none
 
      character*(*) aname, pname
      integer length, type, status
 
      pointer (aptr, array)
      integer array(*)
 
      integer len
 
      call mmblklen (aname, pname, aptr, len, status)
      if (status .ne. 0) then
        call mmgetblk (aname, pname, aptr, length, type, status)
      else if (length .gt. len) then
        call mmnewlen (aname, pname, aptr, length, status)
      else
        status = 0
      end if
 
      return
      end
 
 
      subroutine mm_ovall (aname, pname, aptr, need, extra, length,
     &                     type, status)
 
      implicit none
 
      character*(*) aname, pname
      integer need, extra, length, type, status
 
      pointer (aptr, array)
      integer array(*)
 
      if (need .gt. length) then
        if (length .eq. 0) then
          length = need + extra
          call mmgetblk (aname, pname, aptr, length, type, status)
        else
          length = need + extra
          call mmnewlen (aname, pname, aptr, length, status)
        end if
      end if
 
      return
      end
 
 
      subroutine mmcomprs (iholesize, status)
 
      implicit none
 
      integer iholesize, status
 
      status = 0
 
      return
      end
 
!!
!!  2. PUBLIC INQUIRY PROCEDURES
!!
 
      subroutine mmblklen (aname, pname, aptr, length, status)
 
      implicit none
 
      character*(*) aname, pname
      integer length, status
 
      pointer (aptr, array)
      integer array(*)
 
      call mmfindbk (aname, pname, aptr, length, status)
 
      return
      end
 
 
      subroutine mmfindbk (aname, pname, aptr, length, status)
 
      implicit none
 
      character*(*) aname, pname
      integer length, status
 
      pointer (aptr, array)
      integer array(*)
 
      include 'mm2000.h'
 
      integer index
 
      call mm_index_by_name (aname, pname, index)
      if (index .gt. 0) then
        aptr   = ad_addr(2,index)
        length = ad_addr(5,index)
        status = 0
      else
        aptr = 0
        length = 0
        status = -16
      end if
 
      return
      end
 
 
      subroutine mmgetindex (aname, pname, index, status)
 
      implicit none
 
      character*(*) aname, pname
      integer index, status
 
      call mm_index_by_name (aname, pname, index)
 
      if (index .gt. 0) then
        status = 0
      else
        status = -16
      end if
 
      return
      end
 
 
      subroutine mmgetlen (aptr, length, status)
 
      implicit none
 
      integer length, status
 
      pointer (aptr, array)
      integer array(*)
 
      include 'mm2000.h'
 
      integer index
 
      call mm_index_by_addr (aptr, index)
      if (index .gt. 0) then
        length = ad_addr(5,index)
        status = 0
      else
        status = -16
      end if
 
      return
      end
 
 
      subroutine mmgetnam (aptr, aname, pname, status)
 
      implicit none
 
      character*(*) aname, pname
      integer status
 
      pointer (aptr, array)
      integer array(*)
 
      include 'mm2000.h'
 
      integer index
 
      call mm_index_by_addr (aptr, index)
      if (index .gt. 0) then
        aname = ad_name(1,index)
        pname = ad_name(2,index)
        status = 0
      else
        aname = '$notset'
        pname = '$notset'
        status = -16
      end if
 
      return
      end
 
 
      subroutine mmgetpr (aname, pname, aptr, status)
 
      implicit none
 
      character*(*) aname, pname
      integer status
 
      pointer (aptr, array)
      integer array(*)
 
      include 'mm2000.h'
 
      integer index
 
      call mm_index_by_name (aname, pname, index)
      if (index .gt. 0) then
        aptr   = ad_addr(2,index)
        status = 0
      else
        aptr = 0
        status = -16
      end if
 
      return
      end
 
 
      subroutine mmgettyp (aptr, type, status)
 
      implicit none
 
      integer type, status
 
      pointer (aptr, array)
      integer array(*)
 
      include 'mm2000.h'
 
      integer index
 
      call mm_index_by_addr (aptr, index)
      if (index .gt. 0) then
        type   = ad_addr(4,index)
        status = 0
      else
        type = 0
        status = -16
      end if
 
      return
      end
 
!!
!!  3. PUBLIC DIAGNOSTIC PROCEDURES
!!
 
      subroutine mmverify ()
 
      implicit none
 
      include 'mm2000.h'
 
      pointer (pptr, prolog), (eptr, epilog)
      integer prolog(*), epilog(*)
 
      pointer (aptr, array)
      integer array(*)
 
      integer index, j, length
      character*32 aname, pname
      logical corrupted, stomped_on
 
      index = first
      corrupted = .false.
 
      do while (index .gt. 0)
 
        aname = ad_name(1,index)
        pname = ad_name(2,index)
        pptr  = ad_addr(1,index)
        aptr  = ad_addr(2,index)
        eptr  = ad_addr(3,index)
        length= ad_addr(5,index)
 
        !! CHECK THE PROLOG
        stomped_on = .false.
        do j = 1, AD_SIZE
          if (prolog(j) .ne. ad_addr(j,index)) stomped_on = .true.
        end do
        if (stomped_on) then
          corrupted = .true.
          print *, 'MMVERIFY: Array prolog overwritten'
          print *, '  Array', aname, ' in ', pname
          print *, '  Length=', length, ', Address=', aptr
        end if
 
        !! CHECK THE EPILOG
        stomped_on = .false.
        do j = 1, AD_SIZE
          if (epilog(j) .ne. ad_addr(j,index)) stomped_on = .true.
        end do
        if (stomped_on) then
          corrupted = .true.
          print *, 'MMVERIFY: Array epilog overwritten'
          print *, '  Array', aname, ' in ', pname
          print *, '  Length=', length, ', Address=', aptr
        end if
 
        index = link(index)
      end do
 
      if (corrupted) then
        call mmprint (6)
        call mm_panic ('MMVERIFY', 'Array storage was corrupted')
      end if
 
      return
      end
 
 
      subroutine mmprint (iout)
 
      implicit none
 
      include 'mm2000.h'
 
      integer j, iout
 
      write(iout,8010)
      j = first
      do while (j .gt. 0)
        write(iout,8020) j, ad_addr(5,j), ad_addr(4,j), ad_addr(2,j),
     &                   ad_name(1,j), ad_name(2,j)
        j = link(j)
      end do
 
 8010 format(/,'INDEX',2x,'LENGTH',2x,'TYPE',1x,'ADDRESS',2x,
     &       'NAME',29x,'PARTITION')
 8020 format(i4,1x,i10,1x,i1,' 0x',z8.8,1x,a,1x,a)
 
      return
      end
 
!!
!!  4. PRIVATE PROCEDURES -- DO NOT USE EXTERNALLY
!!
 
      subroutine mm_new_ad (index, status)
 
      implicit none
 
      integer index, status
 
      include 'mm2000.h'
 
      integer max_ad_orig, i, j, mallocf, reallocf
      integer NUM_AD_START, NUM_AD_INCR
      parameter (NUM_AD_START = 500, NUM_AD_INCR = 200)
      external mallocf, reallocf
 
      if (next .eq. 0) then ! No free AD; allocate more.
 
        if (max_ad .gt. 0) then
 
          max_ad_orig = max_ad
          max_ad = max_ad + NUM_AD_INCR
          adnptr = reallocf (adnptr, max_ad*2*BYTES_PER_CHAR)
          adaptr = reallocf (adaptr, max_ad*AD_SIZE*BYTES_PER_PTR)
          lptr   = reallocf (lptr,   max_ad*BYTES_PER_INT)
 
        else
 
          max_ad_orig = 0
          max_ad = NUM_AD_START
          adnptr = mallocf (max_ad*2*BYTES_PER_CHAR)
          adaptr = mallocf (max_ad*AD_SIZE*BYTES_PER_PTR)
          lptr   = mallocf (max_ad*BYTES_PER_INT)
 
        end if
 
        if (adnptr.lt.1 .or. adaptr.lt.1 .or. lptr.lt.1) then
          status = -6
          call mm_panic('MM_NEW_AD',
     &                  'Unable to allocate array desciptors')
        end if
 
        !! Put new AD onto free list
        next = max_ad_orig + 1
        do j = max_ad_orig + 1, max_ad - 1
          link(j) = j + 1
        end do
        link(max_ad) = 0
 
        !! Set new AD to null values
        do j = max_ad_orig + 1, max_ad
          ad_name(1,j) = ' '
          ad_name(2,j) = ' '
          do i = 1, AD_SIZE
            ad_addr(i,j) = 0
          end do
        end do
 
      end if
 
      !! Move AD from free list onto active list
      index = next
      next = link(index)
      link(index) = first
      first = index
      num_ad = num_ad + 1
 
      status = 0
 
      return
      end
 
 
      subroutine mm_index_by_name (aname, pname, index)
 
      implicit none
 
      character*(*) aname, pname
      integer index ,i
 
      include 'mm2000.h'
 
      index = first
      do i=1,100000000
        if (index .lt. 1) then
          index = 0
          return
        end if
        if (aname .eq. ad_name(1,index) .and.
     &      pname .eq. ad_name(2,index)) return
        index = link(index)
      end do
 
      return
      end
 
 
      subroutine mm_index_by_addr (aptr, index)
 
      implicit none
 
      integer index ,i
 
      pointer (aptr, array)
      integer array(*)
 
      include 'mm2000.h'
 
      index = first
      do i=1,100000000
        if (index .lt. 1) then
          index = 0
          return
        end if
        if (aptr .eq. ad_addr(2,index)) return
        index = link(index)
      end do
 
      return
      end
 
 
      subroutine mm_panic (proc, mesg)
      character*(*) proc, mesg
      write(*,'(3a)') proc, ': PANIC! ', mesg
      stop
      end
 
 
      subroutine mm_warn (proc, mesg)
      character*(*) proc, mesg
      write(*,'(3a)') proc, ': Warning! ', mesg
      return
      end
