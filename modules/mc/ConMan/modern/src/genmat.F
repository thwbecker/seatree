      subroutine genmat ()
c
c----------------------------------------------------------------------
c
c  This routine inputs the element material properties.
c
c----------------------------------------------------------------------
c
      use materials
      
      implicit double precision (a-h,o-z)
c
#include<common.h>
c
c.... initialization

      do i = 1, numnp
      dens(i) = zero ; diff(i)   = zero ; tref(i) = zero
      enddo
      do i = 1, numel
       mat(i) = 0
      enddo
      do i = 1, numat
       visc(i) = zero ; alam(i) = zero
       ra(i)   = zero ; dmhu(i) = zero ; tcon(1,i) = zero
       tcon(2,i) = zero
      enddo

      write (iout,1001)  ititle  

c.... read viscosity (non-dimensionalized value)
c
1     read (iin,*,err=1,end=999)  (visc(i),i=1,numat)
      write (iout,1000)  (visc(i),i=1,numat)
c
c.... read penalty number (lamda, non-dimensonalized value)
c
2     read (iin,*,err=2,end=999)  (alam(i),i=1,numat)
      write (iout,1100) (alam(i),i=1,numat)
c      
c.... heat diffusivity (non-dimensionalized value)
c         
3     read (iin,*,err=3,end=999) (diff(i) , i = 1 , numat)
      write (iout,1200)  (diff(i),i=1,numat)       
c
c.... rayleigh number (for initial values for all elements)
c     
5     read (iin,*,err=5,end=999) ( ra(i) , i = 1 , numat)
      write (iout,1300)  (ra(i),i=1,numat)      
c
c.... internal heating paramenter .
c
6     read (iin,*,err=6,end=999) ( dmhu(i) , i = 1 , numat)
      write (iout,1400)  (dmhu(i),i=1,numat)      
c
c.... viscosity law constant
c
7     read (iin,*,err=7,end=999) ( tcon(1,i) , i = 1 , numat)
      write (iout,1500)  (tcon(1,i),i=1,numat)      
c
c.... activation energy for viscosity law
c
8     read (iin,*,err=8,end=999) ( tcon(2,i) , i = 1 , numat)
      write (iout,1600)  (tcon(2,i),i=1,numat)
c
c.... end of file error handling
c
      return
999   call error ('genmat  ','end file',iin)
c
1001  format(' ',80a1,//,
     & ' M a t e r i a l____P r o p e r t i e s                ',1p,//)
1000  format(1p,/
     & ' viscosity_(non_dimensionalized)...........(visc__)=',6e15.8//)
1100  format(1p,/
     & ' lambda_(penalty_number)...................(lamda_)=',6e15.8//)
1200  format(1p,/
     & ' heat_conductivity_(non_dimensionalized)...(diff__)=',6e15.8//)
1300  format(1p,/
     & ' rayleigh_number...........................(ra____)=',6e15.8//)
1400  format(1p,/
     & ' internal_heating_parameter_(J/Kg.s).......(dmhu__)=',6e15.8//)
1500  format(1p,/
     & ' activation_energy_for_viscosity_(J/mol)...(tcon1_)=',6e15.8//)
1600  format(1p,/
     & ' activation_volume_for_viscosity_(m3/mol)..(tcon2_)=',6e15.8//)
c
      end
