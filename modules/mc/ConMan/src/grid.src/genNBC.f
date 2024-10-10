      subroutine genNBC(nb)
c
c-----------------------------------------------------------------------
c subroutine to read type a and b nodes for nuusselt number smoother
c
c        |    |    |    |
c        |    |    |    | 
c   -->  .----.----.----. <-- type b nodes
c        |    |    |    |
c        |    |    |    |
c   -->  .----.----.----. <-- type a nodes
c             bottom
c
c                         top
c             -->  .----.----.----. <-- type a nodes
c                  |    |    |    |
c                  |    |    |    | 
c             -->  .----.----.----. <-- type b nodes
c                  |    |    |    |
c                  |    |    |    |
c 
c-----------------------------------------------------------------------
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
      include '../common.h'
c
      dimension nb(2,1)
c
c.... type a nodes
c
      ncount=0
100   continue
1     read(iin,*,err=1,end=999) nsta,neda,inca
      if (nsta .eq. 0) go to 205
      do 200 i=nsta,neda,inca
        ncount=ncount+1
        nb(1,ncount)=i
200   continue
c
      go to 100
c
c.... type b nodes
c
205   continue
      ncount=0
210   continue
c
2     read(iin,*,err=2,end=999) nstb,nedb,incb
      if (nstb .eq. 0) go to 305
      do 300 i=nstb,nedb,incb
        ncount=ncount+1
        nb(2,ncount)=i
300   continue
c
      go to 210
c
305   continue
      if (necho .eq. 1) then
        write(iout,1000)
        do 400 i=1,nodebn
          write(iout,1001) i,nb(1,i),nb(2,i)
400     continue
      end if
c
      return
c
c.... end of file error handling
c
999   call error('genNBC  ','end file',iin)
1000  format(' smooth boundary nodes ',//,10x,'  #  ',10x,
     &       'nodea',10x,'nodeb')
1001  format(3(10x,i5))
2000  format(' smooth corner nodes ',//,10x,'  #  ',10x,
     &       'nodea',10x,'nodeb',10x,'nodec',10x,'noded')
2001  format(5(10x,i5))
      end
