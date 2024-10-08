ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
c
c  Function to prompt user for an integer value. If user hits return 
c    'getint' is set to 'dfault'.
c
ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      function getint (string, dfault)
c
      character*(*) string
      character*20 tmp

      integer dfault, value, getint

      value = dfault

      write (6, "(a,' [',i8,'] ',$)") string, value
      read (5, "(a20)") tmp
      if (tmp(1:1) .ne. ' ') read (tmp, *) value

      getint = value

      return
      end
