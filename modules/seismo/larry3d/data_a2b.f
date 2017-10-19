c
c       convert data from ascii to binary and do some coordinate 
c	checks
c
c
	character*80 phase	
	real*4 dtime,e1,e2,e3,s1,s2
	integer k,j
	read*,phase
	do k=1,80
	   if(phase(k:k).eq." ")goto11
	enddo
11	open(921,file='data.'//phase(1:k-1)//'.bin',
     &	access='direct',recl=4,status='new',form='unformatted')
	open(422,file='quakes.'//phase(1:k-1)//'.bin',
     &	access='direct',recl=12,status="new",form='unformatted')
	open(522,file='receiv.'//phase(1:k-1)//'.bin',
     &	access='direct',recl=8,status="new",form='unformatted')
	open(20,file="data."//phase(1:k-1)//".ascii",status="old")
	open(21,file="quakes."//phase(1:k-1)//".ascii",status="old")
	open(22,file="receiv."//phase(1:k-1)//".ascii",status="old")
	j=1
	k=1
1	read(20,*,end=2)dtime
	read(21,*)e1,e2,e3
	if(e1.lt.-180)e1=e1+360.
	if(e1.gt.180)e1=e1-360.
	read(22,*)s1,s2
	if(s1.lt.-180)s1=s1+360.
	if(s1.gt.180)s1=s1-360.
	if((s2.gt.-90).and.(s2.lt.90))then
c	no polar receivers
	   write(921,rec=j)dtime
	   write(422,rec=j)e1,e2,e3
	   write(522,rec=j)s1,s2
	   j=j+1
	endif
	k=k+1
	goto1
2	print*,"read ",j-1,"wrote ", k-1, " data"
	close(20)
	close(21)
	close(22)
	close(921)
	close(522)
	close(422)
	end
