      program bin2ascii
      implicit none
      !VARIABLES
      !Parameters
      integer*4,parameter:: nrec = 1000000000       !Largest possible
                                                    !number of records
      !Floating point variables
      real*8,dimension(:),allocatable:: buf8        !Double precision array 
                                                    !to read and write
      real*4,dimension(:),allocatable:: buf4        !Single precision array 
                                                    !to read and write
      !Integer variables
      integer*4,dimension(:),allocatable:: ibuf     !Integer array to 
                                                    !read and write
      integer*4:: ibin                              !Handle for conversion 
                                                    !direction
      integer*4:: reclen                            !Record length
      integer*4:: ncol                              !Number of columns
      !Character variables
      character*132:: filin                         !Input file
      character*132:: filout                        !Output file
      character*1:: type                            !Flag for variable 
                                                    !type
      !Auxiliary variables
      integer*4:: irec                              !Record index
      
      !INPUT SESSION
  101 continue
      write(*,*) 'Binary to ASCII (1) or the reverse (2)?'
      read(*,*) ibin
      if (ibin.ne.1.and.ibin.ne.2) then
        write(*,*) 'Flag not accepted, please insert a valid one'
        go to 101
      end if!(ibin.ne.1.and.ibin.ne.2)
      write(*,*) 'Variable type = ?'
      read(*,*) type
      if (type.ne.'r'.and.type.ne.'i') then
        write(*,*) 'Flag not accepted, please insert a valid one'
      end if!(type.ne.'r'.and.type.ne.'f')
      write(*,*) 'Record length for binary file = ?'
      read(*,*) reclen
      write(*,*) 'Number of columns of ASCII file = ?'
      read(*,*) ncol
      write(*,*) 'Binary file = ?'
      read(*,*) filin
      write(*,*) 'ASCII file = ?'
      read(*,*) filout
      
      !MAIN SESSION
      if (type.eq.'r') then
        if (reclen/ncol.eq.8) then 
          allocate(buf8(ncol))
          open(12,file=filin,access='direct',recl=reclen,&
               form='unformatted')
          open(13,file=filout)
          if (ibin.eq.1) then
            do irec=1,nrec
              read(12,rec=irec,err=901) buf8
              write(13,*) buf8
            end do!irec=1,nrec
          else
            do irec=1,nrec
              read(13,*,end=901) buf8
              write(12,rec=irec) buf8
            end do!irec=1,nrec
          end if!(ibin.eq.1)
        end if!(reclen.eq.8)
        if (reclen/ncol.eq.4) then 
          allocate(buf4(ncol))
            open(12,file=filin,access='direct',recl=reclen,&
               form='unformatted')
          open(13,file=filout)
          if (ibin.eq.1) then
            do irec=1,nrec
              read(12,rec=irec,err=901) buf4
              write(13,*) buf4
            end do!irec=1,nrec
          else
            do irec=1,nrec
              read(13,*,end=901) buf4
              write(12,rec=irec) buf4
            end do!irec=1,nrec
          end if!(ibin.eq.1)
        end if!(reclen.eq.8)
      else
        allocate(ibuf(ncol))
        open(12,file=filin,access='direct',recl=reclen,&
             form='unformatted')
        open(13,file=filout)
        if (ibin.eq.1) then
          do irec=1,nrec
            read(12,rec=irec,err=901) ibuf
            write(13,*) ibuf
          end do!irec=1,nrec
        else
          do irec=1,nrec
            read(13,*,end=901) ibuf
            write(12,rec=irec) ibuf
          end do!irec=1,nrec
        end if!(ibin.eq.1)
      end if!(type.eq.'r')
  901 continue
      end program bin2ascii
