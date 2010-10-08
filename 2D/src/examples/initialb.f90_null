! Initial condition for the vector potential.
! This file contains the expression used for the initial 
! vector potential. You can use temporary real arrays R1-R2 
! of size (n,jsta:jend), and temporary complex arrays C1-C5 
! of size (n,ista:iend) to do intermediate computations. 
! The variable a0 should control the global amplitude of the 
! initial condition, and variables aparam0-9 can be used to 
! control the amplitudes of individual terms. At the end, the 
! potential in spectral space should be stored in the array 
! az (plus bz for the z-component of the magnetic field in 
! 2.5D solvers).

! Null vector potential (2D)
 
      DO j = jsta,jend
         DO i = 1,n
            az(i,j) = 0.0_GP
         END DO
      END DO
