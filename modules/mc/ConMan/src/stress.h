      common /strcom/ ipnb , ipstrs , ippmas

      pointer (ipnb,nb)
      dimension  nb(2,nodebn)

      pointer (ipstrs,stress)
      dimension  stress(6,numnp)

      pointer (ippmas,pmass)
      dimension  pmass(numnp)
