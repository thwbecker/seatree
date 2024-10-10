      character*40  pointer
      character      filfmt*40,         fname*40
c
      pointer = 'ipx'
      ilen = index (pointer, ' ') - 1

      totalsize=0.0
      do num=1,999
      n100  = num / 100
      n100r = mod(num,100)
      n10   = n100r / 10
      n10r  = mod(n100r,10)
      n1    = n10r
      write(filfmt,1000) ilen
1000  format('(a',i3.3,',i1,i1,i1)')
      write(fname,filfmt) pointer(1:ilen),n100,n10,n1
      write(6,*) fname
c  mmgetblk is basically a really fancy wrapper around a malloc call
      write(6,*) "before mmgetblk", num
      call mmgetblk(fname ,'mm', ipx   , 4000000 , 1 , icode)
      write(6,*) "past mmgetblk", num
      size = 8*1*4000000
      totalsize = totalsize+size
      call mmverify()
      enddo
      call mmprint(6)
      write(6,*)totalsize/1.0e6,"MB in", num, "4000000 element arrays"
c
c
      end

