      function kblnk (string)
c
c     returns the number of non-blank characters in string before a blank character
c     is found
c
      character*(*) string
      integer k, kblnk
c
      k = index(string,' ')
      if (k .eq. 0) then
        k = len(string)
      else
        k = k - 1
      endif
c                            
      kblnk = k
      return
      end                    
c
