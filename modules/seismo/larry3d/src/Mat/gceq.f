	subroutine gceq(eplo,epla,stlo,stla,raydel,npts,rayx,rayy,ndim)
c--this routines finds points between epicenter and station for given deltas
	implicit double precision (a-h,o-z)
	dimension urott(3,3),urot(3,3)
	dimension raydel(ndim),rayx(ndim),rayy(ndim)
	radian=90./asin(1.)
c--find Euler angles and rotation matrix
	call pole(epla,eplo,stla,stlo,xlatp,xlonp,azmp,delta)
	alph=xlonp/radian
	beta=(90.-xlatp)/radian
	gama=(180.-azmp-.5*delta)/radian
	call setrot(alph,beta,gama,urot)
c---the inverse of a rotation matrix equals its transpose
	do i=1,3
	do j=1,3
	   urott(i,j)=urot(j,i)
	enddo
	enddo
	do i=1,npts
	   x=raydel(i)
	   y=0.
c---------------------rotate back to unprimed reference
	   call rotll(y,x,y0,x0,urott)
cTEST
	if(x0.lt.0)then
c	   print*,"x0=",x0
	   x0=x0+360.d0
	endif
	   rayx(i)=x0
	   rayy(i)=y0
	enddo
	return
	end

c-----routines from john woodhouse
      subroutine rotll(x,y,x1,y1,urot)
      implicit double precision (a-h,o-z)
      dimension urot(3,3),v(3),v1(3)
      data radian/57.29578/
      a=x/radian
      b=y/radian
      sa=dsin(a)
      ca=dcos(a)
      sb=dsin(b)
      cb=dcos(b)
      v(1)=ca*cb
      v(2)=ca*sb
      v(3)=sa
      do 10 i=1,3
      v1(i)=0.
      do 10 j=1,3
   10 v1(i)=v1(i)+urot(i,j)*v(j)
      x1=datan2(v1(3),dsqrt(dabs(1.-v1(3)*v1(3))))*radian
      y1=datan2(v1(2),v1(1))*radian
      return
      end
      subroutine setrot(alph,beta,gama,urot)
      implicit double precision (a-h,o-z)
      dimension urot(3,3)
      sal=dsin(alph)
      cal=dcos(alph)
      sbe=dsin(beta)
      cbe=dcos(beta)
      sga=dsin(gama)
      cga=dcos(gama)
      urot(1,1)=cga*cbe*cal-sga*sal
      urot(1,2)=cga*cbe*sal+sga*cal
      urot(1,3)=-cga*sbe
      urot(2,1)=-sga*cbe*cal-cga*sal
      urot(2,2)=-sga*cbe*sal+cga*cal
      urot(2,3)=sga*sbe
      urot(3,1)=sbe*cal
      urot(3,2)=sbe*sal
      urot(3,3)=cbe
      return
      end
c find the pole of the path xlatp, xlonp and the azimuth
c at the pole of the epicentre (I thnk).
      subroutine pole(epla,eplo,stla,stlo,xlap,xlop,azmp,delta)
      implicit double precision (a-h,o-z)
c      real epla,eplo,stla,stlo,xlap,xlop,azmp,delta
      data radian/57.295779513082d0/
      th1=(90.d0-epla)/radian
      th2=(90.d0-stla)/radian
      ph1=eplo/radian
      ph2=stlo/radian
      cth1=dcos(th1)
      cth2=dcos(th2)
      sth1=dsin(th1)
      sth2=dsin(th2)
      cph1=dcos(ph1)
      cph2=dcos(ph2)
      sph1=dsin(ph1)
      sph2=dsin(ph2)
      cph21=cph1*cph2+sph1*sph2
      sph21=sph2*cph1-sph1*cph2
      cdel=sth1*sth2*cph21+cth1*cth2
      sdel=dsqrt(1.d0-cdel**2)
      delta=datan2(sdel,cdel)*radian
      ccth=sth1*sth2*sph21/sdel
      scth=dsqrt(1.d0-ccth**2)
      cth=datan2(scth,ccth)
      coco=1.d0/(sdel*scth)
      scph=coco*(cth1*sth2*cph2-cth2*sth1*cph1)
      ccph=coco*(sth1*cth2*sph1-sth2*cth1*sph2)
      cph=datan2(scph,ccph)
      xlap=90.d0-cth*radian
      xlop=cph*radian
      v1=sth1*cph1+sth2*cph2
      v2=sth1*sph1+sth2*sph2
      v3=cth1+cth2
      vn=dsqrt(v1**2+v2**2+v3**2)
      v1=v1/vn
      v2=v2/vn
      v3=v3/vn
      cthm=v3
      sthm=dsqrt(1.d0-v3**2)
      cphm=v1/sthm
      sphm=v2/sthm
      cphmp=ccth*sthm*(cphm*ccph+sphm*scph)-scth*cthm
      sphmp=sthm*(sphm*ccph-cphm*scph)
      phmp=datan2(sphmp,cphmp)
      azmp=180.d0-phmp*radian
      return
      end
