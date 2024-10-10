ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
c
c  Function to choose velocity boundary conditions for corners.
c
ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      function corner (iva, ivb, ifre, ifix)
c
      integer value, corner
c
      if ((iva .eq. ifix) .or. (ivb .eq. ifix)) then
        value = ifix
      else
        value = ifre
      end if
c
      corner = value
c
      return
      end
