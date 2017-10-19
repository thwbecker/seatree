	subroutine linint(yin,yout,xin,xout,n,npt)
	implicit real*8 (a-h,o-z)
	dimension xout(2*n),xin(n),yout(2*n),yin(n)
	j=0
	do i=1,npt-1
	   j=j+1
	   xout(j)=xin(i)
	   yout(j)=yin(i)
	   j=j+1
	   xout(j)=(xin(i+1)+xin(i))/2.d0
	   yout(j)=(yin(i+1)-yin(i))/(xin(i+1)-xin(i))
	   yout(j)=yout(j)*(xout(j)-xin(i))+yin(i)
	enddo
	xout(2*npt-1)=xin(npt)
	yout(2*npt-1)=yin(npt)
	return
	end
	subroutine linint_x(yin,yout,xin,xout,n,npt)
	implicit real*8 (a-h,o-z)
	dimension xout(2*n),xin(n),yout(2*n),yin(n)
	j=0
	do i=1,npt-1
	   j=j+1
	   xout(j)=xin(i)
	   yout(j)=yin(i)
	   if(((360.-yin(i)).lt.90.).and.((360.-yin(i+1)).gt.270.))then
	      y0=yin(i)-360.
	   elseif(((360.-yin(i)).gt.270.).and.((360.-yin(i+1)).lt.90.))then
	      y0=yin(i)+360.
	   else
	      y0=yin(i)
	   endif
	   j=j+1
	   xout(j)=(xin(i+1)+xin(i))/2.d0
	   yout(j)=(yin(i+1)-y0)/(xin(i+1)-xin(i))
	   yout(j)=yout(j)*(xout(j)-xin(i))+y0
	enddo
	xout(2*npt-1)=xin(npt)
	yout(2*npt-1)=yin(npt)
	return
	end
