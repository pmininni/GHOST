! Initial condition for the passive/active scalar/density.
! The scalar can be passive (e.g., in the PHD solver) or 
! active (as in the Boussinesq or compressible solvers).
! This file contains the expression used for the initial 
! concentration of the scalar. You can use temporary 
! real arrays R1-R3 of size (1:nx,1:ny,ksta:kend) and temporary 
! complex arrays C1-C8 of size (1:nz,1:ny,ista:iend) to do 
! intermediate computations. The variable c0 should control 
! the global amplitude of the concentration, and variables 
! cparam0-9 can be used to control the amplitudes of individual 
! terms. At the end, the initial concentration of the scalar in 
! spectral space should be stored in the array th.

! Constant uniform scalar/density
!     c0 : initial uniform scalar/density concentration

! Set scalar/density to zero
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               th(k,j,i) = 0.0_GP
            END DO
         END DO
      END DO

! Add a constant scalar/density
      IF (myrank.eq.0) THEN         
         th(1,1,1) = c0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
         th(1,1,1) = c0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
         th(1,1,1) = c0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
      ENDIF
