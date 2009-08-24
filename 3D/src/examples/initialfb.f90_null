! External electromotive forcing.
! This file contains the expression used for the external 
! electromotive forcing. You can use temporary real arrays 
! R1-R3 of size (1:n,1:n,ksta:kend) and temporary complex 
! arrays C1-C8 of size (n,n,ista:iend) to do intermediate 
! computations. The variable m0 should control the global 
! amplitude of the forcing, and variables mparam0-9 can be
! used to control the amplitudes of individual terms. At the
! end, the three components of the forcing in spectral
! space should be stored in the arrays mx, my, and mz.

! Null electromotive forcing

      DO i = ista,iend
         DO j = 1,n
            DO k = 1,n
               mx(k,j,i) = 0.
               my(k,j,i) = 0.
               mz(k,j,i) = 0.
            END DO
         END DO
      END DO
