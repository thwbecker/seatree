      subroutine span(nsq,xlamin,xlamax,xlomin,xlomax,
     &	nsqrs,nsqtot,nlatzones,eq_incr)
c--finds coordinate range of square number nsq
      implicit real*8 (a-h,o-z)
      dimension nsqrs(nlatzones),nsqtot(nlatzones+1)
      lazone=2
      do while (nsq.gt.nsqtot(lazone))
         lazone=lazone+1
      enddo
      lazone=lazone-1
      nnsq=nsq-nsqtot(lazone)
      xlamin=90.-lazone*eq_incr
      xlamax=xlamin+eq_incr
      grsize=360./nsqrs(lazone)
      xlomax=nnsq*grsize
      xlomin=xlomax-grsize
      return
      end
