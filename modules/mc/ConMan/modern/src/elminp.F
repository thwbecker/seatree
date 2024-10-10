      subroutine elminp ()
c
c----------------------------------------------------------------------
c
c  This routine inputs the element data 
c
c----------------------------------------------------------------------
c
c
      use solver
      use materials
      use variables
      use fe
      use bc
      use mstres

      implicit double precision (a-h,o-z)
c
      include 'common.h'
c
c.... read the material properties
c
      allocate (dens(numnp),  diff(numnp), tref(numnp),
     &         mat(numel), visc(numat), alam(numat),
     &          ra(numat), dmhu(numat), tcon(2,numat), STAT = IERR)
      if (IERR .ne. 0) then
        write(*,*) 'allocation of materials properties fails'
        STOP
      endif
     
      ntotal = ntotal + numat * 7 * 2 + numel
      call genmat ()
        
c
c.... read and set up IEN array
c
      allocate (ien(numel,nen), lmv(numel,nen*ndof), lmt(numel,nen),
     &          STAT = IERR)
      if (IERR .ne. 0) then
        write(*,*) 'allocation of ien, lmv and lmt fails'
        STOP      
      endif

      do i = 1, numel
        do j = 1, nen
          ien(i, j) = 0 ; lmt(i, j) = 0
        enddo
      enddo
      do i = 1, numel
        do j = 1, nen*ndof
          lmv(i,j) = 0
        enddo
      enddo

      ntotal = ntotal + numel*nen*(ndof+2)  
c
c.... diag arrays are for the solver
c
      allocate (idiag(neqv), STAT = IERR)
      if (IERR .ne. 0) then
        write(*,*) 'allocation of idiag fails'
        STOP      
      endif

      do i = 1, neqv
        idiag(i) = 0
      enddo

      ntotal = ntotal + neqv

      if ( (isolve .eq. 2) .or. (isolve .eq. 3) ) then   ! Implicit and Picard
        allocate (idiagt(neqt), STAT = IERR)
        if (IERR .ne. 0) then
          write(*,*) 'allocation of neqt fails'
          STOP      
        endif
        do i = 1, neqt
          idiagt(i) = 0
        enddo
        ntotal = ntotal + neqt
        call genienIP (ien, idiag , mat, lmv , lmt , idv, 
     &                 idiagt  , idt )
      elseif (isolve .eq. 1) then                       ! Explicit
        call genien (ien, idiag , mat, lmv , lmt , idv)
      endif     
c
c.... generate the element shape functions and 
c.... calculate min element dimension per node
c
      allocate (shdx(numel,nen,nipt), shdy(numel,nen,nipt), 
     &          det(numel,nipt), shl(nen,nipt), STAT = IERR)
      if (IERR .ne. 0) then
        write(*,*) 'allocation of shdx, shdy, det and shl fails'
        STOP      
      endif
      
      do i = 1, nen
        do j = 1, nipt
          do iel = 1, numel
            shdx(iel,i,j) = zero ; shdy(iel,i,j) = zero 
            det(iel,j) = zero 
          enddo
          shl(i,j) = zero
        enddo
      enddo
      
      ntotal = ntotal + nen * nipt * 4 * 2 
c
c solver memory
c     
      allocate (vrhs(neqv), vbcr(neqv), trhs(numnp), STAT = IERR)
      if (IERR .ne. 0) then
        write(*,*) 'allocation of vrhs, vbcr and trhs fails'
        STOP      
      endif
      do i = 1, neqv
        vrhs(i) = zero ; vbcr(i) = zero
      enddo
      do i = 1, numnp
        trhs(i) = zero
      enddo
      ntotal = ntotal + 2*neqv*2 + numnp*2
      if ( (isolve .eq. 2) .or. (isolve .eq. 3) ) then   ! Implicit and Picard
        allocate (tlhsa(nsizet), tlhsb(nsizet), STAT = IERR)
        if (IERR .ne. 0) then
          write(*,*) 'allocation of tlhsa and tlhsb fails'
          STOP      
        endif
        do i = 1, nsizet
          tlhsa(i) = zero ; tlhsb(i) = zero
        enddo
        ntotal = ntotal + 2*nsizet*2
      else if (isolve .eq. 1) then                     ! Explicit
        allocate (tlhs(numnp), STAT = IERR)
        if (IERR .ne. 0) then
          write(*,*) 'allocation of tlhs fails'
          STOP      
        endif
        do i = 1, numnp
          tlhs(i) = zero
        enddo
        ntotal = ntotal + numnp*2
      endif 
c arrays for dynamic stress
      allocate (pmass(numnp), stress(5,numnp), STAT=IERR) 
      if (IERR .ne. 0) then
        write(*,*) 'allocation of pmass and stress fails'
        STOP      
      endif
      do i = 1, numnp
        pmass(i) = zero
        do j = 1, 5
          stress(j,i) = zero
        enddo
      enddo
      ntotal = ntotal + 5*numnp*2
c
c allocate memory for solvers
c
      allocate ( vlhsa(nsize), vlhsb(nsize), STAT = IERR)
      if (IERR .ne. 0) then
        write(*,*) 'allocation of vlhs for unfact/unvback fails'
        STOP      
      endif
      do i = 1, nsize
        vlhsa(i) = zero
        vlhsb(i) = zero
      enddo
      ntotal = ntotal + nsize*2
c
c.... return
c
      return
      end
