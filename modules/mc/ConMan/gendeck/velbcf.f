ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
c
c  Function to prompt user for a velocity boundary condition. If user 
c    hits return 'velbcf' is set to 'dfault'.
c
ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      function velbcf (string, dfault, ifre, ifix)
c
      character*(*) string
      character*200 tmp_str
c
      integer velbcf, dfault, value
c
      logical yes
c
      if (dfault .eq. ifre) then
         if (yes (string (1:kblnk (string)) // ' = Unconstrained', 'y')) 
     &    then
          value = ifre
        else
          value = ifix
        end if
      else
        if (yes (string (1:kblnk (string)) // ' = Pinned       ', 'y')) 
     &    then
          value = ifix
        else
          value = ifre
        end if
      end if
c
      velbcf = value
c
      return
      end
