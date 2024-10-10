      common /solcom/ ipvlhs , ipvrhs , ipidiag, ipvbcr,
     &                iptlhsa, iptlhsb, iptrhs,  ipidiat, iptlhs

      pointer (ipvlhs,vlhs)
      dimension  vlhs(nsize)

      pointer (ipvrhs,vrhs)
      dimension  vrhs(neqv)

      pointer (ipvbcr,vbcr)
      dimension  vbcr(neqv)

      pointer (ipidiag,idiag)
      dimension  idiag(neqv)

#ifdef IMPLICIT
      pointer (iptlhsa,tlhsa)
      dimension  tlhsa(nsizet)

      pointer (iptlhsb,tlhsb)
      dimension  tlhsb(nsizet)

      pointer (ipidiat,idiagt)
      dimension  idiagt(neqt)
#else
      pointer (iptlhs,tlhs)
      dimension  tlhs(nEGnp)
#endif

      pointer (iptrhs,trhs)
      dimension  trhs(nEGnp)
