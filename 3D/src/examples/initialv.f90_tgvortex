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

! Superposition of Taylor-Green vortices
!     kdn : minimum wave number
!     kup : maximum wave number

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,n
            DO k = 1,n
               vz(k,j,i) = 0.0_GP
            END DO
         END DO
      END DO

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,n
            DO i = 1,n

            R1(i,j,k) = 0.0_GP
            R2(i,j,k) = 0.0_GP

            DO ki = INT(kdn),INT(kup)
               R1(i,j,k) = R1(i,j,k)+SIN(2*pi*ki*(real(i,kind=GP)-1)/ &
                          real(n,kind=GP))*COS(2*pi*ki*(real(j,kind=GP)-1)/  &
                          real(n,kind=GP))*COS(2*pi*ki*(real(k,kind=GP)-1)/  &
                          real(n,kind=GP))
               R2(i,j,k) = R2(i,j,k)-COS(2*pi*ki*(real(i,kind=GP)-1)/ &
                          real(n,kind=GP))*SIN(2*pi*ki*(real(j,kind=GP)-1)/  &
                          real(n,kind=GP))*COS(2*pi*ki*(real(k,kind=GP)-1)/  &
                          real(n,kind=GP))
            END DO

            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,R1,vx,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R2,vy,MPI_COMM_WORLD)
! We do not normalize if the want the velocity to be given 
! just by sine, cosine with amplitude 1
!      CALL normalize(vx,vy,vz,u0,1,MPI_COMM_WORLD)
