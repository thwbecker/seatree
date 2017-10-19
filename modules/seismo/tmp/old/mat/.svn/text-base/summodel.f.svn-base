       subroutine summodel0(delta,phi,theta,val)
c
c      delta --- epicentral distance ( in degree)
c      phi,theta  --- CAS location (in radian)
c
	parameter (lmax=4,kmax=2)
	dimension a(0:kmax,0:lmax,0:lmax),b(0:kmax,0:lmax,0:lmax)
	dimension f(0:kmax),cmphi(0:lmax),smphi(0:lmax)
	dimension ad(0:lmax,0:lmax),bd(0:lmax,0:lmax)
       dimension p(0:lmax,0:lmax),lmask(0:lmax)
       data init /1/
       save a,b,lmask
       
       
	if(init.eq.1)then
	   lmask(0)=1
	   lmask(1)=0
	   lmask(2)=1
	   lmask(3)=0
	   lmask(4)=1
	   io = 81
	   	   
	   open(io,file='/home/tempo/su/models/inner-core/aniresmdl.dat')
	   read(io,*)maxdeg             ! skip one line

	   if(maxdeg.ne.lmax)stop "error in summodel"
	   
	   write(*,"('degree mask=',132i1)")(lmask(i),i=0,maxdeg)

	   do k=0,kmax
	      do l=0,lmax
	         do m=0,lmax
	            a(k,l,m)=0.0
	            b(k,l,m)=0.0
	         enddo
	      enddo
	   enddo
	   
	   do k=0,kmax
	      do l=0,lmax
	        if(lmask(l).ne.0)then
	           read(io,*)a(k,l,0),(a(k,l,m),b(k,l,m),m=1,l)
	        endif
	      enddo
	   enddo
	
	   close(io)
	   write(*,"(30('-'))")
	   write(*,"('  k  l  m        A         B')")
	   do k=0,kmax
	      do l=0,lmax
	        if(lmask(l).ne.0)then
	          write(*,"(3i3,f10.2)")k,l,0,a(k,l,m)
	          do m=1,l
	            write(*,"(3i3,2f10.2)")k,l,m,a(k,l,m),b(k,l,m)
	          enddo
	        endif
	      enddo
	   enddo
	   init = 0
	   lmask(0)=0       ! do not use l=0 term
	   write(*,"(30('-'))")
	endif
c
c.... calculate cos(m*phi) and sin(m*phi)
c
       cp=cos(phi)
       sp=sin(phi)
       cm=1.0
       sm=0.0
       do m=0,lmax
          cmphi(m)=cm
          smphi(m)=sm
          tm=cm*cp-sm*sp
          sm=sm*cp+cm*sp
          cm=tm
	enddo
c
c.....calculate P(l,m)
c
      call ylm(theta,lmax,p,lmax)		
c
c.... sum in radial direction
c
	x=(180.0-delta)/60.0	
       f(0)=1.0
       f(1)=x*x
       f(2)=x*f(1)
       do l=0,lmax
	   do m=0,l
	      ad(l,m)=0.0
	      bd(l,m)=0.0
	      do k=0,kmax
	         ad(l,m)=ad(l,m)+a(k,l,m)*f(k)
	         bd(l,m)=bd(l,m)+b(k,l,m)*f(k)
	      enddo
	   enddo
	enddo
	
	call sumcoeff0(lmax,p,smphi,cmphi,ad,bd,lmask,val)
	end
	
	subroutine sumcoeff0(lmax,p,smphi,cmphi,ad,bd,lmask,val)
       dimension p(0:lmax,0:lmax),lmask(0:lmax)
	dimension cmphi(0:lmax),smphi(0:lmax)
	dimension ad(0:lmax,0:lmax),bd(0:lmax,0:lmax)
	
	val = 0.0
	do l=0,lmax
	   if(lmask(l).ne.0)then
	      val=val+ad(l,0)*p(l,0)
	      do m=1,l
	         val=val+(ad(l,m)*cmphi(m)+bd(l,m)*smphi(m))*p(l,m)
	      enddo
	   endif
	enddo
	
	return
	end
      subroutine cartesian(s,x)
c.......................................................................
c     This subroutine is to convert the spherical coordinates to
c     cartesian coordinates
c.....Input:
c           s(3) (r,theta,phi) -- spherical coordinates, (in radians)
c.....Output:
c           x(3) (x,y,z)       -- resulting cartesian coordinates
c.......................................................................
      dimension s(3),x(3)
      
      r     = s(1)
      theta = s(2)
      phi   = s(3)
      ct=cos(theta)
      st=sin(theta)
      cp=cos(phi)
      sp=sin(phi)
      x(1)=r*st*cp
      x(2)=r*st*sp
      x(3)=r*ct
      return
      end
      subroutine getrans3(s1,s2,t)
c
c to compute the transformation matrix t, which converts the first 
c point to be on  the x axis and the second point to be on the x-y
c plane on the  new coordinate.
c
c    input:  s1(3) --- first point in spherical coordinate (r,theta,phi)
c            s2(3) --- second point
c    output: t(3,3)--- transformation matrix
c
      dimension s1(3),s2(3),t(3,3),x1(3),x2(3),x3(3)
      
      call cartesian(s1,x1)
      call cartesian(s2,x2)
      
      call outprd(x1,x2,x3)
      call outprd(x3,x1,x2)
      
      call normcarvec(x1,3)
      call normcarvec(x2,3)
      call normcarvec(x3,3)
      
      do i=1,3
         t(1,i)=x1(i)
         t(2,i)=x2(i)
         t(3,i)=x3(i)
      enddo
      
      return
      end
      
 
      subroutine icpointab(depth,delta,d1,d2,ierr)
      parameter (mdel=60,mdep=9)      
      dimension delt(mdel),dep(mdep),ep(mdel,mdep),op(mdel,mdep)
      data init/1/
      save delt,dep,ep,op
      
      ierr = 0
      
      if(init.eq.1)then      
         io = 101
         write(6,"('open Inner-Core entry point file')")
         open(io,file='/home/tempo/su/models/inner-core/ic_point.dat',status='old')
         init = 0
      
         do i=1,mdep
            read(io,*)dep(i)
            write(*,"(f10.2)")dep(i)
            do j=1,mdel
               read(io,*)delt(j),ep(j,i),op(j,i)
cc               write(*,"(f4.0,2f8.2)")delt(j),ep(j,i),op(j,i)
            enddo
         enddo
      
         write(*,"('total number of ',i2,' depths in the table')")mdep 
         write(*,"('total number of ',i2,' deltas in the table')")mdel      
         close(io)
      endif
      
      if(delta.lt.delt(1).or.delta.gt.delt(mdel))then
         ierr = 1
         d1 = 0.0
         d2 = 0.0
         return
      endif
      
      if(depth.lt.0.0.or.depth.gt.950.0)then
         ierr = 2
         d1 = 0.0
         d2 = 0.0
         return
      endif
      
      do i=1,mdep-1
         if(depth.gt.dep(i).and.depth.le.dep(i+1))then
            ndep=i
            goto 1
         endif
      enddo
      
      ndep = mdep
      
1     continue 

      do i=1,mdel-1
         if(delta.ge.delt(i).and.delta.le.delt(i+1))then        
            ndel = i
            goto 2
         endif
      enddo
      
2     d1 = ep(ndel,ndep)
      d2 = op(ndel,ndep)
      return
      end
c.......................................................................
c     nvector.f
c     This subroutine will normalize a vector of length n
c..Input:  x(n)
c..Outpue: x(n)
c.......................................................................
      subroutine normcarvec(x,n)
      dimension x(n)
      
      sum=0.0
      do  i=1,n
         sum=sum+x(i)*x(i)
      enddo
      sum=sqrt(sum)
      do i=1,n
         x(i)=x(i)/sum
      enddo
      
      return
      end
      subroutine transform1(t,x10,x2)
c
c     new coordinate under a transformation
c
      dimension t(3,3),x1(3),x2(3),x10(3)
      
      x1(1)=x10(1)
      x1(2)=x10(2)
      x1(3)=x10(3)

      x2(1)=t(1,1)*x1(1)+t(1,2)*x1(2)+t(1,3)*x1(3)
      x2(2)=t(2,1)*x1(1)+t(2,2)*x1(2)+t(2,3)*x1(3)
      x2(3)=t(3,1)*x1(1)+t(3,2)*x1(2)+t(3,3)*x1(3)
      
      return
      end
      subroutine transform2(t,x10,x2)
c
c     new coordinate under a transformation
c
      dimension t(3,3),x1(3),x2(3),x10(3)
      
      x1(1)=x10(1)
      x1(2)=x10(2)
      x1(3)=x10(3)
      
      x2(1)=t(1,1)*x1(1)+t(2,1)*x1(2)+t(3,1)*x1(3)
      x2(2)=t(1,2)*x1(1)+t(2,2)*x1(2)+t(3,2)*x1(3)
      x2(3)=t(1,3)*x1(1)+t(2,3)*x1(2)+t(3,3)*x1(3)
      
      return
      end
c.......................................................................
c.....ylm.f
c     This subroutine is to compute the spherical harmonics
c     see "Numberical Recipes in C", p194 for the definition of the
c     spherical harmonics.
c..Input:  theta-- spherical coordinate
c          lmax -- angular order
c..Output: p(0:lmax,0:lmax) -- spherical coordinate
c..Note: in this subroutine exp(im*phi) is omitted
c.......................................................................
      subroutine ylm(theta,lmax,p,ld)
      dimension p(0:ld,0:ld)
      data fpi /12.56637061435916/
      
      if(ld.lt.lmax)stop "error in ylm, max dimension exceeded"  
      x=cos(theta)

      y0    =sqrt(1.0/fpi)
      p(0,0)=y0
      p(1,0)=sqrt(3.0/fpi)*x
      fac   =2.0
      sqrx2 =sqrt(1.0-x*x)
      
      do  l=1,lmax
         y0=-sqrt((fac+1.0)/fac)*sqrx2*y0
         fac=fac+2.0
         p(l,l)=y0
         if(l.ne.lmax)p(l+1,l)=sqrt(fac+1.0)*x*y0
      enddo

      do l=2,lmax
         do m=0,l-2
            xl2=2*l
            xlmm=l-m
            xlpm=l+m
            fac =sqrt((xl2+1.0d0)/(xlmm*xlpm))
            fac1=sqrt (xl2-1.0d0)
            fac2=sqrt((xlmm-1)*(xlpm-1)/(xl2-3.0))
            p(l,m)=(x*fac1*p(l-1,m)-fac2*p(l-2,m))*fac
         enddo
      enddo
      return
      end
c.......................................................................
c     This subroutine will compute the outproduct of two vectors (3-d)
c..Input:  x1(3),x2(3)
c..Output: x3(3)
c.......................................................................
      subroutine outprd(x1,x2,x3)
      dimension x1(3),x2(3),x3(3)
      
      x3(1)=x1(2)*x2(3)-x1(3)*x2(2)
      x3(2)=x1(3)*x2(1)-x1(1)*x2(3)
      x3(3)=x1(1)*x2(2)-x1(2)*x2(1)
      return
      end
c.......................................................................
c     this subroutine translates cartesian coordinates to spherical
c     coordinates
c.....input:
c              x, y, z       -- the input cartesian coordinates
c.....output:
c              r, theta, phi -- the resulting spherical coordinates
c                 theta and phi are in radian
c.....note:    when theta is less than 1.0e-6, phi is set to 0.0
c              when r is less than 1.0e-6, phi, theta are set to 0.0
c.......................................................................
      subroutine sphere(x,y,z,r,theta,phi)
      
      r=sqrt(x*x+y*y+z*z)
      theta=0.
      phi=0.
      if(r.le.1.0e-6)return
      theta=acos(z/r)
      if(abs(theta).le.1.0e-6)return
      phi=atan2(y,x)
      if(phi.lt.0.0)phi=phi+2.*3.14159265
      return
      end
