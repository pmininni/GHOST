! Initial condition for the density
! This file contains the expression used for the initial 
! velocity field. You can use temporary real arrays R1-R3 
! of size (1:nx,1:ny,ksta:kend) and temporary complex arrays
! C1-C8 of size (1:nz,1:ny,ista:iend) to do intermediate
! computations. The variable u0 should control the global 
! amplitude of the velocity, and variables vparam0-9 can be
! used to control the amplitudes of individual terms. At the
! end, the three components of the velocity in spectral
! space should be stored in the arrays vx, vy, and vz.

! Set to unity
!     kdn : minimum wave number (rounded to next integer)
!     kup : maximum wave number (rounded to next integer)

!!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
!      DO k = ksta,kend
!!$omp parallel do if (kend-ksta.lt.nth) private (i)
!         DO j = 1,ny
!            DO i = 1,nx
!
!              R1(i,j,k) = rho0
!
!            END DO
!         END DO
!      END DO
!
!      CALL fftp3d_real_to_complex(planrc,R1,rho,MPI_COMM_WORLD)

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               rho(k,j,i) = 0.0_GP
            END DO
         END DO
      END DO

      IF (myrank.eq.0) THEN         
        rho(1,1,1) = rho0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
      ENDIF

