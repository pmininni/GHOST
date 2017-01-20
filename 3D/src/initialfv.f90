! External mechanical forcing.
! This file contains the expression used for the external
! mechanical forcing. You can use temporary real arrays R1-R3
! of size (1:nx,1:ny,ksta:kend) and temporary complex arrays
! C1-C8 of size (1:nz,1:ny,ista:iend) to do intermediate
! computations. The variable f0 should control the global
! amplitude of the forcing, and variables fparam0-9 can be
! used to control the amplitudes of individual terms. At the
! end, the three components of the forcing in spectral
! space should be stored in the arrays fx, fy, and fz.

! Null mechanical forcing

      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               fx(k,j,i) = 0.
               fy(k,j,i) = 0.
               fz(k,j,i) = 0.
            END DO
         END DO
      END DO
