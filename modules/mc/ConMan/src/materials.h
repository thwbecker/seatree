      common /matcom/ ipvisc, ipdiff, ipalam, ipra, ipmhu, iptcon,
     &                ipmat 

      pointer (ipmat,mat)
      dimension  mat(numel)

      pointer (ipvisc,visc)
      dimension  visc(numat)

      pointer (ipdiff,diff)
      dimension  diff(numat)

      pointer (ipra,ra)
      dimension  ra(numat)

      pointer (ipalam,alam)
      dimension  alam(numat)

      pointer (ipmhu,dmhu)
      dimension  dmhu(numat)

      pointer (iptcon,tcon)
      dimension  tcon(3,numat)
