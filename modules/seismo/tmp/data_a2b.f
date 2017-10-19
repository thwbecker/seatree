	character*80 phase	
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
1	read(20,*,end=2)dtime
	read(21,*)e1,e2,e3
	read(22,*)s1,s2
	write(921,rec=j)dtime
	write(422,rec=j)e1,e2,e3
	write(522,rec=j)s1,s2
	j=j+1
	goto1
2	print*,"I have read ",j-1," data"
	close(20)
	close(21)
	close(22)
	close(921)
	close(522)
	close(422)
	end
