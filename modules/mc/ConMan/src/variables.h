      common /varcom/ ipx   , ipv   , ipt   , iptdot, ipdum

      pointer (ipx,x)
      dimension x(nsd,numnp)

      pointer (ipv,v)
      dimension v(ndof,numnp)

      pointer (ipt,t)
      dimension t(numnp)

      pointer (iptdot,tdot)
      dimension tdot(numnp)

      pointer (ipdum,dum)
      dimension dum(numnp)

