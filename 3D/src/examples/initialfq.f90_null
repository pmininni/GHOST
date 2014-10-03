! Thermal forcing for quantum solvers.
! This file contains the expression used for the external
! thermal forcing. You can use temporary real arrays
! R1-R3 of size (1:n,1:n,ksta:kend) and temporary complex
! arrays C1-C8 of size (n,n,ista:iend) to do intermediate
! computations. The variable kttherm can be used to control 
! parameters of the forcing. At the end, the three 
! components of the forcing in spectral space should be 
! stored in the arrays fre and fim.

! Null forcing

      DO i = ista,iend
         DO j = 1,n
            DO k = 1,n
               fre(k,j,i) = 0.
               fim(k,j,i) = 0.
            END DO
         END DO
      END DO
