! External mechanical forcing.
! This file contains the expression used for the external
! mechanical forcing. You can use temporary real arrays
! R1-R3 of size (1:n,1:n,ksta:kend) and temporary complex
! arrays C1-C8 of size (n,n,ista:iend) to do intermediate
! computations. The variable f0 should control the global
! amplitude of the forcing, and variables fparam0-9 can be
! used to control the amplitudes of individual terms. At the
! end, the three components of the forcing in spectral
! space should be stored in the arrays fx, fy, and fz.

! Superposition of Taylor-Green vortices
!     kdn : minimum wave number
!     kup : maximum wave number

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,n
            DO k = 1,n
               fz(k,j,i) = 0.0_GP
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

            DO ki = kdn,kup
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

      CALL fftp3d_real_to_complex(planrc,R1,fx,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R2,fy,MPI_COMM_WORLD)
      CALL normalize(fx,fy,fz,f0,1,MPI_COMM_WORLD)
