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

      DO i = ista,iend
         DO j = 1,n
            DO k = 1,n
               fz(k,j,i) = 0.
            END DO
         END DO
      END DO

      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n

            R1(i,j,k) = 0.
            R2(i,j,k) = 0.

            DO ki = kdn,kup
               R1(i,j,k) = R1(i,j,k)+SIN(2*pi*ki*(float(i)-1)/ &
                          float(n))*COS(2*pi*ki*(float(j)-1)/  &
                          float(n))*COS(2*pi*ki*(float(k)-1)/  &
                          float(n))
               R2(i,j,k) = R2(i,j,k)-COS(2*pi*ki*(float(i)-1)/ &
                          float(n))*SIN(2*pi*ki*(float(j)-1)/  &
                          float(n))*COS(2*pi*ki*(float(k)-1)/  &
                          float(n))
            END DO

            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,R1,fx,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R2,fy,MPI_COMM_WORLD)
      CALL normalize(fx,fy,fz,f0,1,MPI_COMM_WORLD)
