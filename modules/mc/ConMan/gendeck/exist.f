c23456789112345678921234567893123456789412345678951234567896123456789712
c  Steven S. Shapiro
c  1 Nov. 1991
c  4 Sept. 1994
c  Subroutine to determine whether "filename" exists.
c23456789112345678921234567893123456789412345678951234567896123456789712

      function lexist (filename)

      logical lexist
      character*(*) filename

      if (nblen(filename) .eq. 0) then
         lexist = .false.
      else
         inquire (file = filename, exist = lexist)
         if (.not. lexist) then
            print*, filename (1:kblnk (filename)), ' does not exist.'
         else
         end if
      end if

      return

      end
