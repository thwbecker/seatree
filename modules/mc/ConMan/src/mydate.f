      subroutine mydate()
c
c  This appends the date to title,  a useful thing to have in 
c  your output.  However everybody has their own date command.  
c  This keeps the mucking around all in one place.
c
c  SDK 6/18/90
c
      include 'common.h'
c
c GENERIC version
c GENERIC         character*n cdate
c GENERIC         equivalence (cdate, ititle(80-n))
c GENERIC         call date_routine(cdate)
c CRAY    version
c CRAY            character*8 cdate
c CRAY            equivalence (cdate, ititle(72))
c CRAY            call date(cdate)
c SUN     version
c     character*24 cdate
c     equivalence (cdate, ititle(56))
c     call fdate(cdate)
c
      return
      end

