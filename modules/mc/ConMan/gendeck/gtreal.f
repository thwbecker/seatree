ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
c
c  Function to prompt user for a real value. If user hits return 'gtreal' 
c    is set to 'dfault'.
c
ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      function gtreal (string, dfault)
c
      character*(*) string
      character*20 tmp
c
      value = dfault
c
      write (6, "(a,' [', 1pe10.3,'] ',$)") string, value
      read (5, "(a20)") tmp
      if (tmp(1:1) .ne. ' ') read (tmp, *) value
c
      gtreal = value
c
      return
      end
