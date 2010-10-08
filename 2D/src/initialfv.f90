! External mechanical forcing.
! This file contains the expression used for the external 
! mechanical forcing (in streamfunction or velocity formulation 
! depending on the solver). You can use the temporary real array 
! R1 of size (n,jsta:jend), and temporary complex arrays C1, C2 
! of size (n,ista:iend) to do intermediate computations. The 
! variable f0 should control the global amplitude of the forcing, 
! and variables fparam0-9 can be used to control the amplitudes 
! of individual terms. At the end, the forcing in spectral
! space should be stored in the array fk (for the streamfunction, 
! with fz for the z-component in 2.5D solvers), or its components 
! in arrays fx, fy, fz (for velocity field solvers).

! Null forcing in streamfunction formulation (2D)

      DO i = ista,iend
         DO j = 1,n
            fk(j,i) = 0.0_GP
         END DO
      END DO
