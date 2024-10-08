ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
c
c... Function return user's choice of 'yes' or 'no'. Function returns
c       yes = .true. for reposnses which are not understood.
c
ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      function yes (string, default)
c
      logical yes
      character*(*) string
      character*1 default, tmp

      write (*, "(a,' [',a1,']? ',$)") string, default
      read (*, "(a)") tmp
      if (tmp (1:1) .eq. ' ') then
         if ((default .eq. 'n') .or. (default .eq. 'N')) then
            yes = .false.
         else
            yes = .true.
         end if
      else if ((tmp (1:1) .eq. 'n') .or. (tmp (1:1) .eq. 'N')) then
         yes = .false.
      else
         yes = .true.
      end if
c
      return
      end
