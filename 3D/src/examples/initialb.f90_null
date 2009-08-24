! Initial condition for the vector potential.
! This file contains the expression used for the initial 
! vector potential. You can use temporary real arrays R1-R3 
! of size (1:n,1:n,ksta:kend) and temporary complex arrays 
! C1-C8 of size (n,n,ista:iend) to do intermediate 
! computations. The variable a0 should control the global 
! amplitude of the forcing, and variables aparam0-9 can be  
! used to control the amplitudes of individual terms. At the 
! end, the three components of the potential in spectral 
! space should be stored in the arrays ax, ay, and az.

! Null vector potential

      DO i = ista,iend
         DO j = 1,n
            DO k = 1,n
               ax(k,j,i) = 0.
               ay(k,j,i) = 0.
               az(k,j,i) = 0.
            END DO
         END DO
      END DO
