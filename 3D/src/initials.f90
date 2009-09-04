! Initial condition for the passive scalar.
! This file contains the expression used for the initial 
! concentration of the passive scalar. You can use temporary 
! real arrays R1-R3 of size (1:n,1:n,ksta:kend) and temporary 
! complex arrays C1-C8 of size (n,n,ista:iend) to do 
! intermediate computations. The variable c0 should control 
! the global amplitude of the concentration, and variables 
! cparam0-9 can be used to control the amplitudes of individual 
! terms. At the end, the initial passive scalar in spectral 
! space should be stored in the array th.

! Null pasive scalar

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,n
            DO k = 1,n
               th(k,j,i) = 0.0_GP
            END DO
         END DO
      END DO
