      subroutine genfl1 (temp, ninc, inc, node, numgp, ndest, dest)
c
c----------------------------------------------------------------------
c   This program generates floating-point nodal data via 
c isoparametric interpolation used in genfl routine.
c
c         iopt = 1, generation along a line
c              = 2, generation over a surface
c              = 3, generation over a volume
c
c----------------------------------------------------------------------
c
c
      implicit double precision (a-h,o-z)
c
      include 'common.h'
c
      dimension temp(12,*), ninc(*), inc(*), dest(ndest,*), sh(20)
c
c.... initialization
      do i = 1, 20
        sh(i) = zero
      enddo

      iopt = 3
      if (ninc(3) .eq. 0) iopt = 2
      if (ninc(2) .eq. 0) iopt = 1
c
      dr = zero
      ds = zero
      dt = zero
c
      if (ninc(1) .ne. 0) dr = two / dble(ninc(1))
      if (ninc(2) .ne. 0) ds = two / dble(ninc(2))
      if (ninc(3) .ne. 0) dt = two / dble(ninc(3))
c
      ii = ninc(1) + 1
      jj = ninc(2) + 1
      kk = ninc(3) + 1
c
      ni = node
      nj = node
      nk = node
c
      t = -one
      do 500 k=1,kk
      s = -one
      do 400 j=1,jj
      r = -one
      do 300 i=1,ii
c
      call genfl2 (r, s, t, sh, numgp, iopt)
c
      do 200 ides = 1, ndest
        dest(ides,ni) = zero
        do 100 ig = 1, numgp
          dest(ides,ni) = dest(ides,ni) + temp(ides,ig) * sh(ig)
100     continue
200   continue
c
      ni = ni + inc(1)
      r  = r  + dr
300   continue
c
      nj = nj + inc(2)
      ni = nj
      s  = s  + ds
400   continue
c
      nk = nk + inc(3)
      ni = nk
      t  = t + dt
500   continue
c
c.... return
c
      return
      end
      subroutine genfl2 (r, s, t, sh, n, iopt)
c
c----------------------------------------------------------------------
c  This program computes shape functions for isoparametric generation
c used in genfl1 routine.
c
c----------------------------------------------------------------------
c
c
      implicit double precision (a-h,o-z)
c
      include 'common.h'
      dimension sh(*)
c
c.... jump to the right option
c
      goto (100,200,300) iopt
c
c.... ----> 1D generation <-----
c
100   sh(2) = pt5 * r
      sh(1) = pt5 - sh(2)
      sh(2) = pt5 + sh(2)
c
      if (n .eq. 2) return
c
      sh(3) = one   - r * r
      sh(1) = sh(1) - sh(3) * pt5
      sh(2) = sh(2) - sh(3) * pt5
c
      return
c
c.... ----> 2D generation <-----
c
200   r2    = pt5 * r
      r1    = pt5 - r2
      r2    = pt5 + r2
      s2    = pt5 * s
      s1    = pt5 - s2
      s2    = pt5 + s2
      sh(1) = r1  * s1
      sh(2) = r2  * s1
      sh(3) = r2  * s2
      sh(4) = r1  * s2
c
      if (n .eq. 4) return
c
      r3    = one - r * r
      s3    = one - s * s
      sh(5) = r3 * s1
      sh(6) = s3 * r2
      sh(7) = r3 * s2
      sh(8) = s3 * r1
      sh(1) = sh(1) - pt5 * (sh(5) + sh(8))
      sh(2) = sh(2) - pt5 * (sh(6) + sh(5))
      sh(3) = sh(3) - pt5 * (sh(7) + sh(6))
      sh(4) = sh(4) - pt5 * (sh(8) + sh(7))
c
      return
c
c
c.... ----> 3D generation <-----
c
300   r2 = pt5 * r
      r1 = pt5 - r2
      r2 = pt5 + r2
      s2 = pt5 * s 
      s1 = pt5 - s2
      s2 = pt5 + s2
      t2 = pt5 * t
      t1 = pt5 - t2
      t2 = pt5 + t2
c
      rs1   = r1 * s1
      rs2   = r2 * s1
      rs3   = r2 * s2
      rs4   = r1 * s2
      sh(1) = t1 * rs1 
      sh(2) = t1 * rs2
      sh(3) = t1 * rs3
      sh(4) = t1 * rs4
      sh(5) = t2 * rs1
      sh(6) = t2 * rs2
      sh(7) = t2 * rs3
      sh(8) = t2 * rs4
c
      if (n .eq. 8) return
c
      r3     = one - r * r
      s3     = one - s * s
      t3     = one - t * t
c
      sh(17) = t3 * rs1
      sh(18) = t3 * rs2
      sh(19) = t3 * rs3
      sh(20) = t3 * rs4
c
      rs1    = r3 * s1
      rs2    = s3 * r2
      rs3    = r3 * s2
      rs4    = s3 * r1
c
      sh( 9) = t1 * rs1
      sh(10) = t1 * rs2
      sh(11) = t1 * rs3
      sh(12) = t1 * rs4
      sh(13) = t2 * rs1
      sh(14) = t2 * rs2
      sh(15) = t2 * rs3
      sh(16) = t2 * rs4
c
      sh(1)  = sh(1) - pt5 * (sh( 9) + sh(12) + sh(17))
      sh(2)  = sh(2) - pt5 * (sh( 9) + sh(10) + sh(18))
      sh(3)  = sh(3) - pt5 * (sh(10) + sh(11) + sh(19))
      sh(4)  = sh(4) - pt5 * (sh(11) + sh(12) + sh(20))
      sh(5)  = sh(5) - pt5 * (sh(13) + sh(16) + sh(17))
      sh(6)  = sh(6) - pt5 * (sh(13) + sh(14) + sh(18))
      sh(7)  = sh(7) - pt5 * (sh(14) + sh(15) + sh(19))
      sh(8)  = sh(8) - pt5 * (sh(15) + sh(16) + sh(20))
c
      return
c
c.... end of routine
c
      end
