      common /bccom/ ipidv , ipidt , ipvbc , iptbc , ipisuf , ipsuf

      pointer (ipidv,idv)
      dimension  idv(ndof,numnp)

      pointer (ipvbc,vbc)
      dimension  vbc(ndof,numnp)

      pointer (ipidt,idt)
      dimension  idt(numnp)

      pointer (iptbc,tbc)
      dimension  tbc(numnp)

      pointer (ipisuf,isuf)
      dimension  isuf(2,numsuf)

      pointer (ipsuf,suf)
      dimension  suf(3,numsuf)
