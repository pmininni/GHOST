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

! Initialize for QG-balanced initial conditions: set fz=0 (hydrostatic),
! randomize the 2D components, fx, fy s.t. kx fx + ky fy = 0.

!     kdn : minimum wave number (rounded to next integer)
!     kup : maximum wave number (rounded to next integer)

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               fz(k,j,i) = 0.0_GP
            END DO
         END DO
      END DO

      IF (ista.eq.1) THEN
        DO j=2,ny/2+1
           fx     (1,j,1) = 0.0_GP
           fx(1,ny-j+2,1) = 0.0_GP
           fy     (1,j,1) = 0.0_GP
           fy(1,ny-j+2,1) = 0.0_GP
           IF (ky(j).gt.tiny.and.(kk2(1,j,1).le.kup**2).and.(kk2(1,j,1).ge.kdn**2)) THEN
              dump = 1.0_GP  !1./sqrt(kk2(1,j,1))**5
              phase = 2*pi*randu(seed)
              fx     (1,j,1) = (COS(phase)+im*SIN(phase))*dump
              fx(1,ny-j+2,1) = conjg(fx(1,j,1))
              fy     (1,j,1) = -kx(1)/ky(j) * fx(1,j,1)
              fy(1,ny-j+2,1) = conjg(fy(1,j,1))
              
           ENDIF
        ENDDO

        DO k=2,nz/2+1
           fx     (k,1,1) = 0.0_GP
           fx(nz-k+2,1,1) = 0.0_GP
           fy     (k,1,1) = 0.0_GP
           fy(nz-k+2,1,1) = 0.0_GP
           IF (ky(1).gt.tiny.and.(kk2(k,1,1).le.kup**2).and.(kk2(k,1,1).ge.kdn**2)) THEN
              dump = 1.0_GP  !1./sqrt(kk2(k,1,1))**5
              phase = 2*pi*randu(seed)
              fx     (k,1,1) = (COS(phase)+im*SIN(phase))*dump
              fx(nz-k+2,1,1) = conjg(fx(k,1,1))
              fy     (k,1,1) = -kx(1)/ky(1) * fx(k,1,1)
              fy(nz-k+2,1,1) = conjg(fy(k,1,1))
           ENDIF
        ENDDO

        DO j = 2,ny
           DO k = 2,nz/2+1
              fx          (k,j,1) = 0.0_GP
              fx(nz-k+2,ny-j+2,1) = 0.0_GP
              fy          (k,j,1) = 0.0_GP
              fy(nz-k+2,ny-j+2,1) = 0.0_GP
              IF ((kk2(k,j,1).le.kup**2).and.(kk2(k,j,1).ge.kdn**2)) THEN
                  dump = 1.0_GP !1./sqrt(kk2(k,j,1))**5
                  phase = 2*pi*randu(seed)
                  fx        (k,j,1) = (COS(phase)+im*SIN(phase))*dump
                  fx(nz-k+2,ny-j+2,1) = conjg(fx(k,j,1))
                  fy        (k,j,1) = -kx(1)/ky(j) * fx(k,j,1)
                  fy(nz-k+2,ny-j+2,1) = conjg(fy(k,j,1))
              ENDIF
           ENDDO
        ENDDO

        DO i = 2,iend
           DO j = 1,ny
              DO k = 1,nz
                fx(k,j,i) = 0.0_GP
                fy(k,j,i) = 0.0_GP
                IF (ky(j).gt.tiny.and.(kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  dump = 1.0_GP  !1./sqrt(kk2(k,j,i))**5
                  phase = 2*pi*randu(seed)
                  fx(k,j,i) = (COS(phase)+im*SIN(phase))*dump
                  fy(k,j,i) = -kx(i)/ky(j) * fx(k,j,i)
                ENDIF
              ENDDO
           ENDDO
        ENDDO


      ELSE
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
                 fx(k,j,i) = 0.
                 fy(k,j,i) = 0.
                 IF (ky(j).gt.tiny.and.(kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                   dump = 1.0_GP !1./sqrt(kk2(k,j,i))**5
                   phase = 2*pi*randu(seed)
                   fx(k,j,i) = (COS(phase)+im*SIN(phase))*dump
                   fy(k,j,i) = -kx(i)/ky(j) * fx(k,j,i)
                  ENDIF
               ENDDO
            ENDDO
        ENDDO
      ENDIF

      CALL normalize(fx,fy,fz,f0,1,MPI_COMM_WORLD)


