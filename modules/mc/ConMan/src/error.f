      subroutine error (routin, variab, num)
c
c----------------------------------------------------------------------
c
c  This utility routine prints out the error and stops the system
c
c input:
c   routin      : name of the routine where the error occurred
c   variab      : an 8-character error message
c   num         : any integer number associated with the error
c
c----------------------------------------------------------------------
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      include 'common.h'
c
      character*8 routin, variab
c
      data ierchk /0/
c
c.... check for redundant error
c
      if (ierchk .ne. 0) call exit (istat)
      ierchk = 1
c
c.... print the error
c
      write (iout,1000) ititle, routin, variab, num
      if (num .ne. 0) write (iout,1000) ititle, routin, variab, num
      if (num .eq. 0) write (iout,1000) ititle, routin, variab
c
c.... print the time and memory statistics and exit
c
      call timer ('end     ')
      call mmprint(iout)
c
c.... if within processing, print the last time step
c
        if (iflow .eq. 1)
     &    call print(x,v,t,ndof,nelx,nelz,numnp,lstep,time)
c
c.... halt the process
c
      call exit (istat)
c
1000      format(' ',80a1,//,
     &             ' ****** Error occurred in routine <',a8,'>',/,
     &              '  Error code :',a8,:,' : ',i8,//)
      end
