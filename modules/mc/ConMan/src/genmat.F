      subroutine genmat 
c
c----------------------------------------------------------------------
c
c  This routine inputs the element material properties.
c
c----------------------------------------------------------------------
c
c
      implicit double precision (a-h,o-z)
c
c.... deactivate above card(s) for single precision operation
c
#include<common.h>
#include<materials.h>
c
      write (iout,1001)  ititle  
c
c.... read viscosity
c
1     read (iin,*,err=1,end=999)  (visc(i),i=1,numat)
      write (iout,1000)  (visc(i),i=1,numat)
c
c
c.... read lamda
c
2     read (iin,*,err=2,end=999)  (alam(i),i=1,numat)
      write (iout,1100) (alam(i),i=1,numat)
      
c.... heat conductivity
c
3     read (iin,*,err=3,end=999)(diff(i) , i = 1 , numat)
      write (iout,1200)  (diff(i),i=1,numat)
c
c....rayleigh number .
c
4     read (iin,*,err=4,end=999) ( ra(i) , i = 1 , numat)
      write (iout,1300)  (ra(i),i=1,numat)
c
c
c....internal heating paramenter .
c
5     read (iin,*,err=5,end=999) ( dmhu(i) , i = 1 , numat)
      write (iout,1400)  (dmhu(i),i=1,numat)
c
c....viscosity law constant
c
6     read (iin,*,err=6,end=999) ( tcon(1,i) , i = 1 , numat)
      write (iout,1500)  (tcon(1,i),i=1,numat)
c
c....activation energy for viscosity law
c
7     read (iin,*,err=7,end=999) ( tcon(2,i) , i = 1 , numat)
      write (iout,1600)  (tcon(2,i),i=1,numat)
c
8     read (iin,*,err=8,end=999) ( tcon(3,i) , i = 1 , numat)
      write (iout,1700)  (tcon(3,i),i=1,numat)
c
c.... return
c
      return
c
c.... end of file error handling
c
999   call error ('genmat  ','end file',iin)
c
1001  format(' ',80a1,//,
     & ' M a t e r i a l   P r o p e r t i e s                ',1p,//)
1000  format(1p,/
     & ' viscosity . . . . . . . . . . . . . . . . (visc  )=',6e15.5//)
1100  format(1p,/
     & ' lambda  . . . . . . . . . . . . . . . . . (lamda )=',6e15.5//)
1200  format(1p,/
     & ' heat conductivity . . . . . . . . . . . . (diff  )=',6e15.5//)
1300  format(1p,/
     & ' rayleigh number . . . . . . . . . . . . . (ra    )=',6e15.5//)
1400  format(1p,/
     & ' internal heating parameter  . . . . . . . (dmhu  )=',6e15.5//)
1500  format(1p,/
     & ' activation energy for viscosity law . . . (tcon1 )=',6e15.5//)
1600  format(1p,/
     & ' offset for viscosity law  . . . . . . . . (tcon2 )=',6e15.5//)
1700  format(1p,/
     & ' shear Modulusscosity law  . . . . . . . . (tcon3 )=',6e15.5//)
c
      end
