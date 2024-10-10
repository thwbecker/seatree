      subroutine timer (tcode)
c
c----------------------------------------------------------------------
c
c   This routine keeps track of the CPU statistics.
c
c input:
c  tcode      : timer-codes
c                  name of processes to be timed
c                  including 'begin   ', 'end     ','back    '
c                  which are the control input
c
c---------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
c.... remove above card for single-precision operation
c
c
      include 'common.h'
c
      character*8 tcode
c
      character*8 ccode(20)
      real*4      mytime
c
c.... to add new processes to be timed , add character code here
c
        data ccode /'input   ','process ','factor  ','back-slv',
     &              'output  ','pick_dt ','f_vStiff','f_vRes  ',
     &              'f_tRes  ','predict ','begin   ','end     ',
     &              '        ','        ','        ','        ',
     &              '        ','        ','        ','        '/
      data ibg /11/
      data cpu /20*0.e0/
c
c.... get the cpu time
c
      cpunew = mytime()
c
c.... Find the code
c
      do 100 i = 1, ibg + 1
        ii = i
        if (tcode .eq. ccode(i)) goto 200
100      continue
      call error ('timer   ', tcode, 0)
c
200      continue
c
c.... initialize the timer
c
      if (ii .eq. ibg) then 
c
c.... call the macro routine to read the job status
c
        cpuold = mytime()
c
        icode  = ibg
        return
      endif
c
c.... update the status
c
      cpu(icode)    = cpu(icode)    + cpunew - cpuold
      icode         = ii
      cpuold        = mytime()
      if (ii .ne. ibg+1) return
c
c.... print out
c
      ccode(ibg) = 'total   '
      do 300 i= 1, ibg-1
        cpu(ibg)    = cpu(ibg)    + cpu(i)
300   continue
c
      write (iout,1000) ititle
      write (iout,1100) (ccode(i), cpu(i), i=1,ibg)
c
c.... return
c
      return
c
1000      format(' ',80a1,//,
     &  ' E x e c u t i o n   T i m e   S t a t i s t i c s    ',//,
     &  ' name        CPU time  '  )
1100      format(1x,a8,2x,f10.2)
c
      end
