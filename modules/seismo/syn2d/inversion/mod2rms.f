      character*80 nome
      print*,"file with model?"
      read*,nome
      open(1,file=nome,status="old")
      n=0
      rms=0
 1    read(1,*,end=2)k,x
      rms=rms+x*x
      n=n+1
      goto1
 2    close(1)
      rms=rms/n
      rms=sqrt(rms)
      open(1,file="rms.txt",access="append")
      write(1,*)rms
      close(1)
      end
