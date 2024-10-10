      common /fecomm/ ipien , iplmv , iplmt , ipshl ,
     &                ipshdx, ipshdy, ipdet 

      pointer (ipien,ien)
      dimension  ien(numel,nen)

      pointer (iplmv,lmv)
      dimension  lmv(numel,nen*ndof)

      pointer (iplmt,lmt)
      dimension  lmt(numel,nen)

      pointer (ipshl,shl)
      dimension  shl(nen,nipt)

      pointer (ipshdx,shldx)
      dimension  shldx(nen,nipt)

      pointer (ipshdy,shldy)
      dimension  shldy(nen,nipt)

c     pointer (ipdet,det)
c     dimension  det(nipt)
