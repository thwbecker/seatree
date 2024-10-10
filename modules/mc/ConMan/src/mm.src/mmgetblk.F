c $Log:   /pvcs.config/mc3dvv/src/mmgetblk.F_a  $
c
c   Rev 1.0   Mon Jul 14 14:21:46 1997   pvcs
c 
      subroutine mmgetblk(string1,string2,ipoint,num,idum,icode)
      character*(*) string1,string2
      integer*8 ipoint
    
      nbyte = 8*num
      
      ipoint = mallocf(nbyte)

      return
      end
