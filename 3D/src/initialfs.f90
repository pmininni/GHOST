! External source of passive/active scalar.
! The scalar can be passive (e.g., in the PHD solver) or 
! active (as in the Boussinesq solvers). 
! This file contains the expression used for the external 
! forcing of the scalar quantity in many solvers. You can 
! use temporary real arrays R1-R3 of size (1:n,1:n,ksta:kend) 
! and temporary complex arrays C1-C8 of size (n,n,ista:iend) 
! to do intermediate computations. The variable s0 should 
! control the global amplitude of the forcing, and variables 
! mparam0-9 can be used to control the amplitudes of 
! individual terms. At the end, the forcing in spectral
! space should be stored in the array fs.

! Null source

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,n
            DO k = 1,n
               fs(k,j,i) = 0.0_GP
            END DO
         END DO
      END DO
