! Initial condition for the momentum density. Density
! must be set before computing.
!
! This file contains the expression used for the initial 
! velocity field. You can use temporary real arrays R1-R3 
! of size (1:nx,1:ny,ksta:kend) and temporary complex arrays
! C1-C8 of size (1:nz,1:ny,ista:iend) to do intermediate
! computations. The variable u0 should control the global 
! amplitude of the velocity, and variables vparam0-9 can be
! used to control the amplitudes of individual terms. At the
! end, the three components of the velocity in spectral
! space should be stored in the arrays vx, vy, and vz.

! Superposition of ABC vortices
!     kdn : minimum wave number (rounded to next integer)
!     kup : maximum wave number (rounded to next integer)
!     fparam0: A amplitude
!     fparam1: B amplitude
!     fparam2: C amplitude

      IF ( (abs(Lx-Ly).gt.tinyd).or.(abs(Lx-Lz).gt.tinyd) ) THEN
        IF (myrank.eq.0) &
           PRINT *,'ABC initial conditions require Lx=Ly=Lz'
        STOP
      ENDIF

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               sz(k,j,i) = 0.0_GP
               vz(k,j,i) = 0.0_GP
            END DO
         END DO
      END DO

      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx

            R1(i,j,k) = 0.
            R2(i,j,k) = 0.
            R3(i,j,k) = 0.

            DO ki = INT(kdn),INT(kup)
               R1(i,j,k) = R1(i,j,k)+(fparam1*COS(2*pi*ki*(real(j,kind=GP)-1)/ &
                          real(ny,kind=GP))+fparam2*SIN(2*pi*ki*(real(k,kind=GP)-1)/ &
                          real(nz,kind=GP)))/ki**2
               R2(i,j,k) = R2(i,j,k)+(fparam0*SIN(2*pi*ki*(real(i,kind=GP)-1)/ &
                          real(nx,kind=GP))+fparam2*COS(2*pi*ki*(real(k,kind=GP)-1)/ &
                          real(nz,kind=GP)))/ki**2
               R3(i,j,k) = R3(i,j,k)+(fparam0*COS(2*pi*ki*(real(i,kind=GP)-1)/ &
                          real(nx,kind=GP))+fparam1*SIN(2*pi*ki*(real(j,kind=GP)-1)/ &
                          real(ny,kind=GP)))/ki**2
            END DO

            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,R1,fx,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R2,fy,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R3,fz,MPI_COMM_WORLD)
      CALL normalize(fx,fy,fz,f0,1,MPI_COMM_WORLD)

!#if !defined(DENSITY_)
!  #error "DENSITY_ not defined"
!#endif

#if defined(MOM_)
      ! Dealias:
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               C1(k,j,i) = fx (k,j,i)
               C2(k,j,i) = fy (k,j,i)
               C3(k,j,i) = fz (k,j,i)
               C4(k,j,i) = rho(k,j,i)
               IF (kn2(k,j,i).gt.kmax) THEN
                  C1 (k,j,i) = 0.0_GP
                  C2 (k,j,i) = 0.0_GP
                  C3 (k,j,i) = 0.0_GP
                  C4 (k,j,i) = 0.0_GP
               ENDIF
            END DO
         END DO
      END DO

      CALL fftp3d_complex_to_real(plancr,C1 ,R1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,C2 ,R2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,C3 ,R3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,C4 ,R7,MPI_COMM_WORLD)

      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               R1(i,j,k) = R1(i,j,k)*R7(i,j,k) * tmp
               R2(i,j,k) = R2(i,j,k)*R7(i,j,k) * tmp
               R3(i,j,k) = R3(i,j,k)*R7(i,j,k) * tmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,R1,fx,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R2,fy,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R3,fz,MPI_COMM_WORLD)

#endif
