      subroutine ppmprt(ipsol, ngpars, vfile, nplot)
c
c.... program to print temperature data in ppm format
c
c
c-----------------------------------------------------------------------
c
c input:
c  vfile        : file name
c  nplot        : Imagetool frame number
c
c-----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
      parameter (ipsolx = 301, ipsolz = 101)
      include 'common.h'
c
      dimension red(255) , blue(255) , green(255)
      integer ipsol(ipsolx,ipsolz), ngpars(2)
      character*(*)  vfile
      character      filfmt*120,         fname*120

c.... initialization
      do i = 1, 255
         red(i) = zero ; blue(i) = 255 ; green(i) = 255
      enddo

c
c.... setup the file name format
c
      do 10 i = 1, 33
        red(i) = 0.0D+0
        green(i) = 0.0D+0
        blue(i) = 0.5D+0 + 0.5D+0*float(i-1)/32.0D+0
 10   continue
      do 11 i = 34, 66
        red(i) = 0.0D+0
        green(i) = 0.5D+0*float(i-34)/32.0D+0
        blue(i) = 1.0D+0
 11   continue
      do 12 i = 67, 99
        red(i) =  0.0D+0
        green(i) = 0.5D+0 + 0.5D+0*float(i-67)/32.0D+0
        blue(i) = 1.0D+0
 12   continue
      do 13 i = 100, 132
        red(i) =  0.5D+0*float(i-100)/32.0D+0
        green(i) = 1.0D+0 
        blue(i) = 1.0D+0 - 0.5D+0*float(i-100)/32.0D+0
 13   continue
      do 14 i = 133, 165
        red(i) =  0.5D+0  + 0.5D+0*float(i-133)/32.0D+0
        green(i) = 1.0D+0
        blue(i) = 0.5D+0 - 0.5D+0*float(i-133)/32.0D+0
 14   continue
      do 15 i = 165, 197
        red(i) =  1.0D+0 
        green(i) = 1.0D+0 - 0.5D+0*float(i-165)/32.0D+0
        blue(i) = 0.0D+0
 15   continue
      do 16 i = 198, 230
        red(i) =  1.0D+0 
        green(i) = 0.5D+0 - 0.5D+0*float(i-198)/32.0D+0
        blue(i) = 0.0D+0
 16   continue
      do 17 i = 231, 255
        red(i) =  1.0D+0 - 0.5D+0*float(i-231)/24.0D+0
        green(i) = 0.0D+0 
        blue(i) = 0.0D+0
 17   continue
      ilen = index (vfile, ' ') - 1
      n1000 = nplot / 1000
      n1000r= mod(nplot,1000)
      n100  = n1000r / 100
      n100r = mod(n1000r,100)
      n10   = n100r / 10
      n10r  = mod(n100r,10)
      n1    = n10r
      write(filfmt,1000) ilen
1000  format('(a',i4.4,',i1,i1,i1,i1,''.ppm'',1x)')
c
c.... open file
c
      write(fname,filfmt) vfile(1:ilen),n1000,n100,n10,n1
      open (unit=itufout, file=fname, status='unknown', err=999)
c
      write(itufout,"('P3')")
      write(itufout,*) ngpars(1), ngpars(2)
      write(itufout,"('255')")
      do 101 nx=1,ipsolx
        do 101 nz=1,ipsolz
        if(ipsol(nx,nz) .gt. 255) ipsol(nx,nz) = 255 
 101  continue
  
      do 100 nz=ngpars(2),1,-1
       write(itufout,2000)      
     &    ( 
     &      int(255*  red(ipsol(nx,nz))), 
     &      int(255*green(ipsol(nx,nz))), 
     &      int(255* blue(ipsol(nx,nz))),
     &      nx = 1,ngpars(1),1
     &    )
     
      do 100 nx=1,ipsolx
        if(ipsol(nx,nz) .lt. 0) write(6,*) nx,nz,ipsol(nx,nz)
        if(ipsol(nx,nz) .gt. 255) write(6,*) nx,nz,ipsol(nx,nz)
  100 continue

      call flush(itufout)
c
      close (unit=itufout)
      return
c
 2000 format(301(3i5))
  999 call error ("print 2","('Unable to open file; unit ',i3)",itufout)
      end
