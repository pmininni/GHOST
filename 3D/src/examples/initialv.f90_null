! Initial condition for the velocity.
! This file contains the expression used for the initial
! velocity field. You can use temporary real arrays R1-R3
! of size (1:n,1:n,ksta:kend) and temporary complex arrays
! C1-C8 of size (n,n,ista:iend) to do intermediate
! computations. The variable u0 should control the global
! amplitude of the velocity, and variables vparam0-9 can be
! used to control the amplitudes of individual terms. At the
! end, the three components of the velocity in spectral
! space should be stored in the arrays vx, vy, and vz.

! Null velocity field

      DO i = ista,iend
         DO j = 1,n
            DO k = 1,n
               vx(k,j,i) = 0.
               vy(k,j,i) = 0.
               vz(k,j,i) = 0.
            END DO
         END DO
      END DO
