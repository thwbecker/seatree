      subroutine aniso_cor0(slat,slon,selv,elat,elon,edep,resic,ierr)
c.....to caculate the travel-time residual due to the anisotropy of the
c     inner core
c.....Input:   slat (-90->+90) station latitude
c              slon (-180->+180) station longitude
c              selv (km), station elevation (above sea level as +)
c              elat (-90->+90)  earthquake latitude
c              elon (-180->+180) earthquake longitude
c              edep (km), earthquake depth ( always positive)
c.....Output:  resic (second)
c......................................................................
      parameter (rad=3.14159265/180.0,radius=6371.0)
      parameter (delmin=120.0)
      dimension s1(3),s2(3),x1(3),x2(3),x3(3),xe(3),t(3,3)

      s1(1)=(radius+ selv)/radius            ! station spherical coordinate
      s1(2)=(90.0-slat)*rad                  ! theta
      xlon1=slon
      if(xlon1.lt.0.0)xlon1=xlon1+360.0
      s1(3)=xlon1*rad                        ! phi

      s2(1)=(radius-edep)/radius             ! eq  spherical coordinate
      s2(2)=(90.0-elat)*rad
      xlon1=elon
      if(xlon1.lt.0.0)xlon1=xlon1+360.0
      s2(3)=xlon1*rad
         
      call getrans3(s1,s2,t)                 ! get transformation matrix
      call cartesian(s2,x2)
      call transform1(t,x2,x3)               ! eq location in new coor
      angle=atan2(x3(2),x3(1))/rad           ! epicentral distance      
      call icpointab(edep,angle,d1,d2,ierr)  ! enter & exit location
       
      resic  = 0.0
      if(ierr.eq.0.and.angle.ge.delmin.and.angle.le.180.0)then            
         xe(1)=cos(d1*rad)-cos(d2*rad)
         xe(2)=sin(d1*rad)-sin(d2*rad)
         xe(3)=0.0
            
         call normcarvec(xe,3)
         call transform2(t,xe,x1)
            
         call sphere(x1(1),x1(2),x1(3),r1,theta1,phi1)
         call summodel0(angle,phi1,theta1,resic)
cc         if(angle.ge.150.0)then
cc         write(*,"(f8.2,f8.2,f7.2,f5.2)")angle,phi1*180.0/3.1416,
cc     &theta1*180.0/3.1416,resic
cc        else
cc         write(*,*)
cc         endif
cccc      else
cc         write(*,*)
cc         write(*,"(f8.2,f8.2,f7.2,f5.2)")angle,phi1*180.0/3.1416,
cc     &theta1*180.0/3.1416,99.9
      endif
      return
      end
