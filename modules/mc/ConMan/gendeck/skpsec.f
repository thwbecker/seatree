ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
c
c  Subroutine to read and ignore the contents of every line in a file
c     until "nzero" consecutive entries of value zero are found in one
c     line.
c
ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      function skpsec (iunit, nzero)
c
      dimension entry (20)
c
      parameter (zero = 0.0)
c
      finish = 1.0
c
      do 10 while (finish .ne. zero)
        read (iunit, *, end = 10, err = 10) (entry (i), i = 1, nzero)
        j = 1
        tmp = zero
c
        do 20 while ((tmp .eq. zero) .and. (j .ne. nzero))
          tmp = entry (j)
          j = j + 1
  20    continue
c
        finish = tmp
c
  10  continue
c
      skpsec = 0
      return
      end
