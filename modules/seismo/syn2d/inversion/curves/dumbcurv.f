c--from output of dumbder (first and second derivative of a curve)
c--calculate curvature of the curve
	implicit real*4 (a-h,o-z)
      character*80 infile
      parameter(n=4) !x values, f(x), f'(x), f"(x)
      dimension x(20)
      print*,"input file with x, f(x), f'(x), f''(x)?"
      read*,infile
      open(1,file=infile,status="old")
      do l=1,80
         if(infile(l:l).eq." ")goto2
      enddo
 2    open(2,file=infile(1:l-1)//".crv")
 1    read(1,*,end=3)(x(i),i=1,n)
c--curvature k from http://mathworld.wolfram.com/Curvature.html eq. (11)
      den=1.+(x(3)*x(3))
      den=sqrt(den**3)
      curv=x(4)/den
cTEST
c	print*,x(1),curv,den
      write(2,*)(x(i),i=1,n),curv
      goto1
 3    continue
      close(1)
      close(2)
      end
