! External mechanical forcing.
! This file contains the expression used for the external 
! mechanical forcing. You can use temporary real arrays 
! R1-R3 of size (1:nx,1:ny,ksta:kend) and temporary complex 
! arrays C1-C8 of size (nz,ny,ista:iend) to do intermediate 
! computations. The variable f0 should control the global 
! amplitude of the forcing, and variables fparam0-9 can be
! used to control the amplitudes of individual terms. At the
! end, the three components of the forcing in spectral
! space should be stored in the arrays fx, fy, and fz. 

! Superposition of harmonic modes with random phases 
! correlated to give a tunable total relative helicity 
! [after Pouquet & Patterson, JFM 1978]
!     kdn    : minimum wave number !     kup    : maximum wave number
!     fparam0: sin(2*fparam0) is the total relative helicity
!     fparam8, fparam9 are used for temporal storage

      IF (ista.eq.1) THEN
         C1(1,1,1) = 0.
         C2(1,1,1) = 0.
         C3(1,1,1) = 0.
         DO j = 2,ny/2+1 
            IF ((kk2(1,j,1).le.kup**2).and.(kk2(1,j,1).ge.kdn**2)) THEN
               dump = 1./sqrt(kk2(1,j,1))
               phase = 2*pi*randu(seed)
               C1(1,j,1) = (COS(phase)+im*SIN(phase))*dump
               C1(1,ny-j+2,1) = conjg(C1(1,j,1))
               phase = 2*pi*randu(seed)
               C2(1,j,1) = (COS(phase)+im*SIN(phase))*dump
               C2(1,ny-j+2,1) = conjg(C2(1,j,1))
               phase = 2*pi*randu(seed)
               C3(1,j,1) = (COS(phase)+im*SIN(phase))*dump
               C3(1,ny-j+2,1) = conjg(C3(1,j,1))
            ELSE
               C1(1,j,1) = 0.
               C1(1,ny-j+2,1) = 0.
               C2(1,j,1) = 0.
               C2(1,ny-j+2,1) = 0.
               C3(1,j,1) = 0.
               C3(1,ny-j+2,1) = 0.
            ENDIF

         END DO
         DO k = 2,nz/2+1

            IF ((kk2(k,1,1).le.kup**2).and.(kk2(k,1,1).ge.kdn**2)) THEN
               dump = 1./sqrt(kk2(k,1,1))
               phase = 2*pi*randu(seed)
               C1(k,1,1) = (COS(phase)+im*SIN(phase))*dump
               C1(nz-k+2,1,1) = conjg(C1(k,1,1))
               phase = 2*pi*randu(seed)
               C2(k,1,1) = (COS(phase)+im*SIN(phase))*dump
               C2(nz-k+2,1,1) = conjg(C2(k,1,1))
               phase = 2*pi*randu(seed)
               C3(k,1,1) = (COS(phase)+im*SIN(phase))*dump
               C3(nz-k+2,1,1) = conjg(C3(k,1,1))
            ELSE
               C1(k,1,1) = 0.
               C1(nz-k+2,1,1) = 0.
               C2(k,1,1) = 0.
               C2(nz-k+2,1,1) = 0.
               C3(k,1,1) = 0.
               C3(nz-k+2,1,1) = 0.
            ENDIF

         END DO
         DO j = 2,ny
            DO k = 2,nz/2+1
     
            IF ((kk2(k,j,1).le.kup**2).and.(kk2(k,j,1).ge.kdn**2)) THEN
               dump = 1./sqrt(kk2(k,j,1))
               phase = 2*pi*randu(seed)
               C1(k,j,1) = (COS(phase)+im*SIN(phase))*dump
               C1(nz-k+2,ny-j+2,1) = conjg(C1(k,j,1))
               phase = 2*pi*randu(seed)
               C2(k,j,1) = (COS(phase)+im*SIN(phase))*dump
               C2(nz-k+2,ny-j+2,1) = conjg(C2(k,j,1))
               phase = 2*pi*randu(seed)
               C3(k,j,1) = (COS(phase)+im*SIN(phase))*dump
               C3(nz-k+2,ny-j+2,1) = conjg(C3(k,j,1))
            ELSE
               C1(k,j,1) = 0.
               C1(nz-k+2,ny-j+2,1) = 0.
               C2(k,j,1) = 0.
               C2(nz-k+2,ny-j+2,1) = 0.
               C3(k,j,1) = 0.
               C3(nz-k+2,ny-j+2,1) = 0.
            ENDIF

            END DO
         END DO
         DO i = 2,iend
            DO j = 1,ny
               DO k = 1,nz

               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  dump = 1./sqrt(kk2(k,j,i))
                  phase = 2*pi*randu(seed)
                  C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
               ELSE
                  C1(k,j,i) = 0.
                  C2(k,j,i) = 0.
                  C3(k,j,i) = 0.
               ENDIF

               END DO
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz

               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  dump = 1./sqrt(kk2(k,j,i))
                  phase = 2*pi*randu(seed)
                  C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
               ELSE
                  C1(k,j,i) = 0.
                  C2(k,j,i) = 0.
                  C3(k,j,i) = 0.
               ENDIF

               END DO
            END DO
        END DO
      ENDIF

      CALL rotor3(C2,C3,C4,1)
      CALL rotor3(C1,C3,C5,2)
      CALL rotor3(C1,C2,C6,3)
      CALL normalize(C4,C5,C6,f0,1,MPI_COMM_WORLD)

!C4,C5,C6 are the components of the first random vector (say v1)
!Repeating this procedure we create the second random vector (v2)

      IF (ista.eq.1) THEN
         C1(1,1,1) = 0.
         C2(1,1,1) = 0.
         C3(1,1,1) = 0.
         DO j = 2,ny/2+1

            IF ((kk2(1,j,1).le.kup**2).and.(kk2(1,j,1).ge.kdn**2)) THEN
               dump = 1./sqrt(kk2(1,j,1))
               phase = 2*pi*randu(seed)
               C1(1,j,1) = (COS(phase)+im*SIN(phase))*dump
               C1(1,ny-j+2,1) = conjg(C1(1,j,1))
               phase = 2*pi*randu(seed)
               C2(1,j,1) = (COS(phase)+im*SIN(phase))*dump
               C2(1,ny-j+2,1) = conjg(C2(1,j,1))
               phase = 2*pi*randu(seed)
               C3(1,j,1) = (COS(phase)+im*SIN(phase))*dump
               C3(1,ny-j+2,1) = conjg(C3(1,j,1))
            ELSE
               C1(1,j,1) = 0.
               C1(1,ny-j+2,1) = 0.
               C2(1,j,1) = 0.
               C2(1,ny-j+2,1) = 0.
               C3(1,j,1) = 0.
               C3(1,ny-j+2,1) = 0.
            ENDIF

         END DO
         DO k = 2,nz/2+1

            IF ((kk2(k,1,1).le.kup**2).and.(kk2(k,1,1).ge.kdn**2)) THEN
               dump = 1./sqrt(kk2(k,1,1))
               phase = 2*pi*randu(seed)
               C1(k,1,1) = (COS(phase)+im*SIN(phase))*dump
               C1(nz-k+2,1,1) = conjg(C1(k,1,1))
               phase = 2*pi*randu(seed)
               C2(k,1,1) = (COS(phase)+im*SIN(phase))*dump
               C2(nz-k+2,1,1) = conjg(C2(k,1,1))
               phase = 2*pi*randu(seed)
               C3(k,1,1) = (COS(phase)+im*SIN(phase))*dump
               C3(nz-k+2,1,1) = conjg(C3(k,1,1))
            ELSE
               C1(k,1,1) = 0.
               C1(nz-k+2,1,1) = 0.
               C2(k,1,1) = 0.
               C2(nz-k+2,1,1) = 0.
               C3(k,1,1) = 0.
               C3(nz-k+2,1,1) = 0.
            ENDIF

         END DO
         DO j = 2,ny
            DO k = 2,nz/2+1
     
            IF ((kk2(k,j,1).le.kup**2).and.(kk2(k,j,1).ge.kdn**2)) THEN
               dump = 1./sqrt(kk2(k,j,1))
               phase = 2*pi*randu(seed)
               C1(k,j,1) = (COS(phase)+im*SIN(phase))*dump
               C1(nz-k+2,ny-j+2,1) = conjg(C1(k,j,1))
               phase = 2*pi*randu(seed)
               C2(k,j,1) = (COS(phase)+im*SIN(phase))*dump
               C2(nz-k+2,ny-j+2,1) = conjg(C2(k,j,1))
               phase = 2*pi*randu(seed)
               C3(k,j,1) = (COS(phase)+im*SIN(phase))*dump
               C3(nz-k+2,ny-j+2,1) = conjg(C3(k,j,1))
            ELSE
               C1(k,j,1) = 0.
               C1(nz-k+2,ny-j+2,1) = 0.
               C2(k,j,1) = 0.
               C2(nz-k+2,ny-j+2,1) = 0.
               C3(k,j,1) = 0.
               C3(nz-k+2,ny-j+2,1) = 0.
            ENDIF

            END DO
         END DO
         DO i = 2,iend
            DO j = 1,ny
               DO k = 1,nz

               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  dump = 1./sqrt(kk2(k,j,i))
                  phase = 2*pi*randu(seed)
                  C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
               ELSE
                  C1(k,j,i) = 0.
                  C2(k,j,i) = 0.
                  C3(k,j,i) = 0.
               ENDIF

               END DO
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz

               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  dump = 1./sqrt(kk2(k,j,i))
                  phase = 2*pi*randu(seed)
                  C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
               ELSE
                  C1(k,j,i) = 0.
                  C2(k,j,i) = 0.
                  C3(k,j,i) = 0.
               ENDIF

               END DO
            END DO
        END DO
      ENDIF

      CALL rotor3(C2,C3,C7,1)
      CALL rotor3(C1,C3,C8,2)
      CALL rotor3(C1,C2,C3,3)
      CALL normalize(C7,C8,C3,f0,1,MPI_COMM_WORLD)

! So far, v1 = (C4,C5,C6) and v2 = (C7,C8,C3) are two random 
! normalized vectors in Fourier space. Correlating vectors 
! v1 and v2 we create the forcing vectors fx, fy, and fz.

      fparam8=SIN(fparam0)
      fparam9=COS(fparam0)

! fx: components y and z of u2 are calculated in order to get 
! the x component of curl(u2)

      DO i = ista,iend
          DO j = 1,ny
             DO k = 1,nz
                IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  C1(k,j,i) = fparam8*C5(k,j,i)+fparam9*C8(k,j,i)
                  C2(k,j,i) = fparam8*C6(k,j,i)+fparam9*C3(k,j,i)
                ENDIF
             END DO
          END DO
      END DO

      CALL rotor3(C1,C2,C1,1)

      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  dump = 1./sqrt(kk2(k,j,i))
                  fx(k,j,i) = fparam9*C4(k,j,i)+fparam8*C7(k,j,i)+ &
                              C1(k,j,i)*dump
               ENDIF
            END DO
         END DO
      END DO

! fy: components x and z of u2 are calculated in order to get 
! the y component of curl(u2)

      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  C1(k,j,i) = fparam8*C4(k,j,i)+fparam9*C7(k,j,i)
                  C2(k,j,i) = fparam8*C6(k,j,i)+fparam9*C3(k,j,i)
               ENDIF
            END DO
         END DO
      END DO

      CALL rotor3(C1,C2,C1,2)

      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  dump = 1./sqrt(kk2(k,j,i))
                  fy(k,j,i) = fparam9*C5(k,j,i)+fparam8*C8(k,j,i)+ &
                              C1(k,j,i)*dump
               ENDIF
            END DO
         END DO
      END DO

! fz: components x and y of u2 are calculated in order to get
! the z component of curl(u2)

      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  C1(k,j,i) = fparam8*C4(k,j,i)+fparam9*C7(k,j,i)
                  C2(k,j,i) = fparam8*C5(k,j,i)+fparam9*C8(k,j,i)
               ENDIF
            END DO
         END DO
      END DO

      CALL rotor3(C1,C2,C1,3)

      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  dump = 1./sqrt(kk2(k,j,i))
                  fz(k,j,i) = fparam9*C6(k,j,i)+fparam8*C3(k,j,i)+ &
                              C1(k,j,i)*dump
               ENDIF
            END DO
         END DO
      END DO

      CALL normalize(fx,fy,fz,f0,1,MPI_COMM_WORLD)
