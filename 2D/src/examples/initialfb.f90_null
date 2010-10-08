! External electromotive forcing.
! This file contains the expression used for the external 
! electromotive forcing. You can use temporary real arrays 
! R1-R2 of size (n,jsta:jend), and temporary complex arrays 
! C1-C5 of size (n,ista:iend) to do intermediate computations. 
! The variable m0 should control the global amplitude of the 
! forcing, and variables mparam0-9 can be used to control the 
! amplitudes of individual terms. At the end, the forcing in 
! spectral space should be stored in the array mk (plus mz 
! for the z-component of the magnetic field in 2.5D solvers).

! Null electromotive forcing

 
      DO j = jsta,jend
         DO i = 1,n
            mk(i,j) = 0.0_GP
         END DO
      END DO