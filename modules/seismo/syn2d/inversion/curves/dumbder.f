      character nome*80,chj*2,chk*2
      dimension x(20,2),xd(20)
      print*,"input file?"
      read*,nome
      open(1,file=nome,status="old")
      do l=1,80
         if(nome(l:l).eq." ")goto3
      enddo
 3    print*,"how many columns?"
      read*,n
      print*,"which column to derive?"
      read*,j
      print*,"wrt which column?"
      read*,k
      write(chj,"(i2.2)")j
      write(chk,"(i2.2)")k
      open(2,file=nome(1:l-1)//".der."//chj//".wrt."//chk)
      read(1,*,end=2)(x(i,1),i=1,n)
 1    read(1,*,end=2)(x(i,2),i=1,n)
         dy=x(j,2)-x(j,1)
         dx=x(k,2)-x(k,1)
	 do m=1,n
            xd(m)=(x(m,2)+x(m,1))/2.
	    x(m,1)=x(m,2)
	 enddo
         write(2,*)(xd(i),i=1,n),dy/dx
      goto1
 2    continue
      close(1)
      close(2)
      end
