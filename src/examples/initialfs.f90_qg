! External source of passive/active scalar.
! The scalar can be passive (e.g., in the PHD solver) or 
! active (as in the Boussinesq solvers). 
! This file contains the expression used for the external 
! use temporary real arrays R1-R3 of size (1:nx,1:ny,ksta:kend)
! and temporary complex arrays C1-C8 of size (1:nz,1:ny,ista:iend)
! to do intermediate computations. The variable s0 should 
! control the global amplitude of the forcing, and variables 
! mparam0-9 can be used to control the amplitudes of 
! individual terms. At the end, the forcing in spectral
! space should be stored in the array fs.
!
!
! Scalar initialized with geostrophic condition, assuming 
! hydrostatic balance in z (so, must use TG or something
! like initialv.f90_bal).
!     skdn    : minimum wave number
!     skup    : maximum wave number
!
! In wavespace:
!             theta = -i/N k_z P,
! where
!             P = -1/k_perp^2 * f * (curl v)_z,
! and
!             f = 2 Omega
!
! Must ensure that velocity forcing is set before the
! call to this code is made. This velocity, set, e.g., 
! from initialfv.f90_tg, should set vz=0, and still
! be Div-free.
!
!
! First, compute z-vorticity:
      CALL rotor3(fx,fy,C4,1)  ! C4=(curl f)_z
!
! Compute complex constant, i/N:
      cdump = im/bvfreq
!
      IF (ista.eq.1) THEN
        DO j=2,ny/2+1
           fs    (1,j,1) = 0.0_GP
           fs(1,ny-j+2,1) = 0.0_GP
           tmp = kx(1)**2+ky(j)**2
           IF (tmp.gt.tiny) THEN
              C2    (1,j,1) = -2.0*(omegaz/tmp)*C4(1,j,1) ! balanced pressure
              fs    (1,j,1) = -cdump*kz(k)*C2(1,j,1)     ! theta
              fs(1,ny-j+2,1) = conjg(fs(1,j,1))
           ENDIF
        ENDDO

        DO k=2,nz/2+1
           fs     (k,1,1) = 0.0_GP
           fs(nz-k+2,1,1) = 0.0_GP
           tmp = kx(1)**2+ky(1)**2
           IF (tmp.gt.tiny) THEN
              C2    (k,1,1) = -2.0*(omegaz/tmp)*C4(k,1,1) ! balanced pressure
              fs    (k,1,1) = -cdump*kz(k)*C2(k,1,1)     ! theta
              fs(nz-k+2,1,1) = conjg(fs(k,1,1))
           ENDIF
        ENDDO

        DO j = 2,ny
           DO k = 2,nz/2+1
              fs          (k,j,1) = 0.0_GP
              fs(nz-k+2,ny-j+2,1) = 0.0_GP
              tmp = kx(1)**2+ky(j)**2
              IF (tmp.gt.tiny) THEN
                  C2    (k,j,1) = -2.0*(omegaz/tmp)*C4(k,j,1) ! balanced pressure
                  fs    (k,j,1) = -cdump*kz(k)*C2(k,j,1)     ! theta
                  fs(nz-k+2,ny-j+2,1) = conjg(fs(k,j,1))
              ENDIF
           ENDDO
        ENDDO

        DO i = 2,iend
           DO j = 1,ny
              DO k = 1,nz
                fs(k,j,i) = 0.0_GP
                tmp = kx(i)**2+ky(j)**2
                IF (tmp.gt.tiny) THEN
                  C2    (k,j,i) = -2.0*(omegaz/tmp)*C4(k,j,i) ! balanced pressure
                  fs    (k,j,i) = -cdump*kz(k)*C2(k,j,i)     ! theta
                ENDIF
              ENDDO
           ENDDO
        ENDDO


      ELSE
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
                 tmp = kx(i)**2+ky(j)**2
                  C2(k,j,i) = 0.0_GP
                  IF ( tmp.gt.tiny ) then
                    C2(k,j,i) = -2.0*(omegaz/tmp)*C4(k,j,i) ! balanced pressure
                    fs(k,j,i) = -cdump*kz(k)*C2(k,j,i)     ! theta
                  ENDIF
               ENDDO
            ENDDO
        ENDDO
      ENDIF

      CALL variance(fs,tmp,1)
      CALL MPI_BCAST(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               fs(k,j,i) = fs(k,j,i)*s0/sqrt(tmp)
            END DO
         END DO
      END DO

