!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Extra subroutines to compute the energy spectrum and energy 
! transfer functions in the HD, MHD, and Hall-MHD equations 
! in 3D in the anisotropic case (e.g. in the rotating frame 
! or with an imposed external magnetic field). Quantities 
! parallel and perpendicular to the z direction in Fourier 
! space can be computed. You should use the FFTPLANS and 
! MPIVARS modules (see the file 'fftp_mod.f90') in each 
! program that calls any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2005 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.uba.ar 
!=================================================================

!*****************************************************************
      SUBROUTINE specpara(a,b,c,nmb,kin,hel)
!-----------------------------------------------------------------
!
! Computes the energy and helicity power spectrum in 
! the direction parallel to the preferred direction 
! (rotation or uniform magnetic field). As a result, the 
! k-shells are planes with normal (0,0,kz) (kz=0,...,n/2). 
! The output is written to a file by the first node.
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     nmb: the extension used when writting the file
!     kin: =2 skips energy spectrum computation
!          =1 computes the kinetic spectrum
!          =0 computes the magnetic spectrum
!     hel: =1 computes the helicity spectrum
!          =0 skips helicity spectrum computation
!
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1)            :: Ek,Ektot
      COMPLEX, INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX, DIMENSION(n,n,ista:iend)             :: c1,c2,c3
      REAL                :: tmp
      INTEGER, INTENT(IN) :: kin,hel
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.
      END DO
!
! Computes the curl of the field if needed
!
      IF ((kin.eq.0).or.(hel.eq.1)) THEN
         CALL rotor3(b,c,c1,1)
         CALL rotor3(a,c,c2,2)
         CALL rotor3(a,b,c3,3)
      ENDIF
!
! Computes the kinetic energy spectrum
!
      tmp = 1./float(n)**6
      IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn)
            DO j = 1,n
               DO k = 1,n
                  kmn = int(abs(ka(k))+1)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!$omp atomic
                     Ek(kmn) = Ek(kmn)+                            &
                       (abs(a(k,j,1))**2+abs(b(k,j,1))**2+         &
                       abs(c(k,j,1))**2)*tmp
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!$omp atomic
                        Ek(kmn) = Ek(kmn)+                         &
                          2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+    &
                          abs(c(k,j,i))**2)*tmp
                     ENDIF
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!$omp atomic
                        Ek(kmn) = Ek(kmn)+                         &
                          2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+    &
                          abs(c(k,j,i))**2)*tmp
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the magnetic energy spectrum
!
      ELSE IF (kin.eq.0) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn)
            DO j = 1,n
               DO k = 1,n
                  kmn = int(abs(ka(k))+1)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!$omp atomic
                     Ek(kmn) = Ek(kmn)+                            &
                       (abs(c1(k,j,1))**2+abs(c2(k,j,1))**2+       &
                       abs(c3(k,j,1))**2)*tmp
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!$omp atomic
                        Ek(kmn) = Ek(kmn)+                         &
                          2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+  &
                          abs(c3(k,j,i))**2)*tmp
                     ENDIF
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!$omp atomic
                        Ek(kmn) = Ek(kmn)+                         &
                          2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+  &
                          abs(c3(k,j,i))**2)*tmp
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      IF (kin.le.1) THEN
         CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) THEN
            IF (kin.eq.1) THEN
               OPEN(1,file='kspecpara.' // nmb // '.txt')
            ELSE
               OPEN(1,file='mspecpara.' // nmb // '.txt')
            ENDIF
            WRITE(1,20) Ektot
   20       FORMAT( E23.15 ) 
            CLOSE(1)
         ENDIF
      END IF
!
! Computes the helicity spectrum
!
      IF (hel.eq.1) THEN
         DO i = 1,n/2+1
            Ek(i) = 0.
         END DO
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn)
            DO j = 1,n
               DO k = 1,n
                  kmn = int(abs(ka(k))+1)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!$omp atomic
                     Ek(kmn) = Ek(kmn)+                                 &
                       (real(a(k,j,1)*conjg(c1(k,j,1)))+                &
                       real(b(k,j,1)*conjg(c2(k,j,1)))+                 &
                       real(c(k,j,1)*conjg(c3(k,j,1))))*tmp
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!$omp atomic
                        Ek(kmn) = Ek(kmn)+                              &
                          2*(real(a(k,j,i)*conjg(c1(k,j,i)))+           &
                          real(b(k,j,i)*conjg(c2(k,j,i)))+              &
                          real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
                     ENDIF
                 END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!$omp atomic
                        Ek(kmn) = Ek(kmn)+                              &
                          2*(real(a(k,j,i)*conjg(c1(k,j,i)))+           &
                          real(b(k,j,i)*conjg(c2(k,j,i)))+              &
                          real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
         CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         0,MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) THEN
            IF (kin.eq.1) THEN
               OPEN(1,file='khelipara.' // nmb // '.txt')
            ELSE IF (kin.eq.0) THEN
               OPEN(1,file='mhelipara.' // nmb // '.txt')
            ELSE
               OPEN(1,file='ghelipara.' // nmb // '.txt')
            ENDIF
            WRITE(1,30) Ektot
   30       FORMAT( E23.15 )
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE specpara

!*****************************************************************
      SUBROUTINE specperp(a,b,c,nmb,kin,hel)
!-----------------------------------------------------------------
!
! Computes the energy and helicity power spectrum in the 
! direction perpendicular to the preferred direction 
! (rotation or uniform magnetic field). The k-shells are 
! cylindrical surfaces (kperp=1,...,N/2+1). The output is 
! written to a file by the first node.
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     nmb: the extension used when writting the file
!     kin: =2 skips energy spectrum computation
!          =1 computes the kinetic spectrum
!          =0 computes the magnetic spectrum
!     hel: =1 computes the helicity spectrum
!          =0 skips helicity spectrum computation
!
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1)            :: Ek,Ektot
      COMPLEX, INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX, DIMENSION(n,n,ista:iend)             :: c1,c2,c3
      REAL                :: tmp
      INTEGER, INTENT(IN) :: kin,hel
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.
      END DO
!
! Computes the curl of the field if needed
!
      IF ((kin.eq.0).or.(hel.eq.1)) THEN
         CALL rotor3(b,c,c1,1)
         CALL rotor3(a,c,c2,2)
         CALL rotor3(a,b,c3,3)
      ENDIF
!
! Computes the kinetic energy spectrum
!
      tmp = 1./float(n)**6
      IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn)
            DO j = 1,n
               kmn = int(sqrt(ka(1)**2+ka(j)**2)+.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  DO k = 1,n
!$omp atomic
                     Ek(kmn) = Ek(kmn)+                            &
                       (abs(a(k,j,1))**2+abs(b(k,j,1))**2+         &
                       abs(c(k,j,1))**2)*tmp
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
!$omp atomic
                        Ek(kmn) = Ek(kmn)+                         &
                          2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+    &
                          abs(c(k,j,i))**2)*tmp
                     END DO
                  ENDIF
               END DO
            END DO
          ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
!$omp atomic
                        Ek(kmn) = Ek(kmn)+                         &
                          2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+    &
                          abs(c(k,j,i))**2)*tmp
                     END DO
                  ENDIF
               END DO
            END DO
          ENDIF
!
! Computes the magnetic energy spectrum
!
      ELSE IF (kin.eq.0) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn)
            DO j = 1,n
               kmn = int(sqrt(ka(1)**2+ka(j)**2)+.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  DO k = 1,n
!$omp atomic
                     Ek(kmn) = Ek(kmn)+                            &
                       (abs(c1(k,j,1))**2+abs(c2(k,j,1))**2+       &
                       abs(c3(k,j,1))**2)*tmp
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
!$omp atomic
                        Ek(kmn) = Ek(kmn)+                         &
                          2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+  &
                          abs(c3(k,j,i))**2)*tmp
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
!$omp atomic
                        Ek(kmn) = Ek(kmn)+                         &
                          2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+  &
                          abs(c3(k,j,i))**2)*tmp
                     END DO
                  ENDIF
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      IF (kin.le.1) THEN
         CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) THEN
            IF (kin.eq.1) THEN
               OPEN(1,file='kspecperp.' // nmb // '.txt')
            ELSE
               OPEN(1,file='mspecperp.' // nmb // '.txt')
            ENDIF
            WRITE(1,20) Ektot
   20       FORMAT( E23.15 ) 
            CLOSE(1)
         ENDIF
      END IF
!
! Computes the helicity spectrum
!
      IF (hel.eq.1) THEN
         DO i = 1,n/2+1
            Ek(i) = 0.
         END DO
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn)
            DO j = 1,n
               kmn = int(sqrt(ka(1)**2+ka(j)**2)+.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  DO k = 1,n
!$omp atomic
                     Ek(kmn) = Ek(kmn)+                                 &
                       (real(a(k,j,1)*conjg(c1(k,j,1)))+                &
                       real(b(k,j,1)*conjg(c2(k,j,1)))+                 &
                       real(c(k,j,1)*conjg(c3(k,j,1))))*tmp
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
!$omp atomic
                        Ek(kmn) = Ek(kmn)+                              &
                          2*(real(a(k,j,i)*conjg(c1(k,j,i)))+           &
                          real(b(k,j,i)*conjg(c2(k,j,i)))+              &
                          real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
!$omp atomic
                        Ek(kmn) = Ek(kmn)+                              &
                          2*(real(a(k,j,i)*conjg(c1(k,j,i)))+           &
                          real(b(k,j,i)*conjg(c2(k,j,i)))+              &
                          real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
                     END DO
                  ENDIF
               END DO
            END DO
         ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
         CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         0,MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) THEN
            IF (kin.eq.1) THEN
               OPEN(1,file='kheliperp.' // nmb // '.txt')
            ELSE IF (kin.eq.0) THEN
               OPEN(1,file='mheliperp.' // nmb // '.txt')
            ELSE
               OPEN(1,file='gheliperp.' // nmb // '.txt')
            ENDIF
            WRITE(1,30) Ektot
   30       FORMAT( E23.15 )
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE specperp

!*****************************************************************
      SUBROUTINE entpara(a,b,c,d,e,f,nmb,kin)
!-----------------------------------------------------------------
!
! Computes the energy transfer in the direction parallel 
! to the preferred direction (rotation or uniform magnetic 
! field) in 3D Fourier space. The k-shells are planes with 
! normal (0,0,kz) (kz=0,...,n/2). The output is written to 
! a file by the first node.
!
! Parameters
!     a  : field component in the x-direction
!     b  : field component in the y-direction
!     c  : field component in the z-direction
!     d  : nonlinear term in the x-direction
!     e  : nonlinear term in the y-direction
!     f  : nonlinear term in the z-direction
!     nmb: the extension used when writting the file
!     kin: =0 computes the magnetic energy transfer
!          =1 computes the kinetic energy transfer
!          =2 computes the Lorentz force transfer
!
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1)          :: Ek,Ektot
      COMPLEX, INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX, INTENT(IN), DIMENSION(n,n,ista:iend) :: d,e,f
      REAL                :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.
      END DO
!
! Computes the kinetic energy transfer
!
      tmp = 1./float(n)**6
      IF ((kin.eq.1).or.(kin.eq.2)) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn)
            DO j = 1,n
               DO k = 1,n
                  kmn = int(abs(ka(k))+1)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!$omp atomic
                     Ek(kmn) = Ek(kmn)+                             &
                       (real(a(k,j,1)*conjg(d(k,j,1)))+             &
                       real(b(k,j,1)*conjg(e(k,j,1)))+              &
                       real(c(k,j,1)*conjg(f(k,j,1))))*tmp
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!$omp atomic
                        Ek(kmn) = Ek(kmn)+                          &
                          2*(real(a(k,j,i)*conjg(d(k,j,i)))+        &
                          real(b(k,j,i)*conjg(e(k,j,i)))+           &
                          real(c(k,j,i)*conjg(f(k,j,i))))*tmp
                     ENDIF
                 END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!$omp atomic
                        Ek(kmn) = Ek(kmn)+                          &
                          2*(real(a(k,j,i)*conjg(d(k,j,i)))+        &
                          real(b(k,j,i)*conjg(e(k,j,i)))+           &
                          real(c(k,j,i)*conjg(f(k,j,i))))*tmp
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the magnetic energy transfer
!
      ELSE
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn)
            DO j = 1,n
               DO k = 1,n
                  kmn = int(abs(ka(k))+1)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!$omp atomic
                     Ek(kmn) = Ek(kmn)+ka2(k,j,1)*                  &
                       (real(a(k,j,1)*conjg(d(k,j,1)))+             &
                       real(b(k,j,1)*conjg(e(k,j,1)))+              &
                       real(c(k,j,1)*conjg(f(k,j,1))))*tmp
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!$omp atomic
                        Ek(kmn) = Ek(kmn)+2*ka2(k,j,i)*             &
                          (real(a(k,j,i)*conjg(d(k,j,i)))+          &
                          real(b(k,j,i)*conjg(e(k,j,i)))+           &
                          real(c(k,j,i)*conjg(f(k,j,i))))*tmp
                     ENDIF
                 END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!$omp atomic
                        Ek(kmn) = Ek(kmn)+2*ka2(k,j,i)*             &
                          (real(a(k,j,i)*conjg(d(k,j,i)))+          &
                          real(b(k,j,i)*conjg(e(k,j,i)))+           &
                          real(c(k,j,i)*conjg(f(k,j,i))))*tmp
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF (kin.eq.0) THEN
            OPEN(1,file='mtranpara.' // nmb // '.txt')
         ELSEIF (kin.eq.1) THEN
            OPEN(1,file='ktranpara.' // nmb // '.txt')
         ELSE
            OPEN(1,file='jtranpara.' // nmb // '.txt')
         ENDIF
         WRITE(1,40) Ektot
   40    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE entpara

!*****************************************************************
      SUBROUTINE entperp(a,b,c,d,e,f,nmb,kin)
!-----------------------------------------------------------------
!
! Computes the energy transfer in the direction perpendicular
! to the preferred direction (rotation or uniform magnetic 
! field) in 3D Fourier space. The k-shells are cylindrical 
! surfaces (kperp=1,...,N/2+1). The output is written to a 
! file by the first node.
!
! Parameters
!     a  : field component in the x-direction
!     b  : field component in the y-direction
!     c  : field component in the z-direction
!     d  : nonlinear term in the x-direction
!     e  : nonlinear term in the y-direction
!     f  : nonlinear term in the z-direction
!     nmb: the extension used when writting the file
!     kin: =0 computes the magnetic energy transfer
!          =1 computes the kinetic energy transfer
!          =2 computes the Lorentz force transfer
!
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1)            :: Ek,Ektot
      COMPLEX, INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX, INTENT(IN), DIMENSION(n,n,ista:iend) :: d,e,f
      REAL                :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.
      END DO
!
! Computes the kinetic energy transfer
!
      tmp = 1./float(n)**6
      IF ((kin.eq.1).or.(kin.eq.2)) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn)
            DO j = 1,n
               kmn = int(sqrt(ka(1)**2+ka(j)**2)+.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  DO k = 1,n
!$omp atomic
                     Ek(kmn) = Ek(kmn)+                             &
                       (real(a(k,j,1)*conjg(d(k,j,1)))+             &
                       real(b(k,j,1)*conjg(e(k,j,1)))+              &
                       real(c(k,j,1)*conjg(f(k,j,1))))*tmp
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
!$omp atomic
                        Ek(kmn) = Ek(kmn)+                          &
                          2*(real(a(k,j,i)*conjg(d(k,j,i)))+        &
                          real(b(k,j,i)*conjg(e(k,j,i)))+           &
                          real(c(k,j,i)*conjg(f(k,j,i))))*tmp
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
!$omp atomic
                        Ek(kmn) = Ek(kmn)+                          &
                          2*(real(a(k,j,i)*conjg(d(k,j,i)))+        &
                          real(b(k,j,i)*conjg(e(k,j,i)))+           &
                          real(c(k,j,i)*conjg(f(k,j,i))))*tmp
                     END DO
                  ENDIF
               END DO
            END DO
         ENDIF
!
! Computes the magnetic energy transfer
!
      ELSE
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn)
            DO j = 1,n
               kmn = int(sqrt(ka(1)**2+ka(j)**2)+.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  DO k = 1,n
!$omp atomic
                     Ek(kmn) = Ek(kmn)+ka2(k,j,1)*                  &
                       (real(a(k,j,1)*conjg(d(k,j,1)))+             &
                       real(b(k,j,1)*conjg(e(k,j,1)))+              &
                       real(c(k,j,1)*conjg(f(k,j,1))))*tmp
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
!$omp atomic
                        Ek(kmn) = Ek(kmn)+2*ka2(k,j,i)*             &
                          (real(a(k,j,i)*conjg(d(k,j,i)))+          &
                          real(b(k,j,i)*conjg(e(k,j,i)))+           &
                          real(c(k,j,i)*conjg(f(k,j,i))))*tmp
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
!$omp atomic
                        Ek(kmn) = Ek(kmn)+2*ka2(k,j,i)*             &
                          (real(a(k,j,i)*conjg(d(k,j,i)))+          &
                          real(b(k,j,i)*conjg(e(k,j,i)))+           &
                          real(c(k,j,i)*conjg(f(k,j,i))))*tmp
                     END DO
                  ENDIF
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF (kin.eq.0) THEN
            OPEN(1,file='mtranperp.' // nmb // '.txt')
         ELSEIF (kin.eq.1) THEN
            OPEN(1,file='ktranperp.' // nmb // '.txt')
         ELSE
            OPEN(1,file='jtranperp.' // nmb // '.txt')
         ENDIF
         WRITE(1,40) Ektot
   40    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE entperp

!*****************************************************************
      SUBROUTINE heltpara(a,b,c,d,e,f,nmb,kin)
!-----------------------------------------------------------------
!
! Computes the helicity transfer in the direction parallel 
! to the preferred direction (rotation or uniform magnetic 
! field) in 3D Fourier space. The k-shells are planes with 
! normal (0,0,kz) (kz=0,...,n/2). The output is written to 
! a file by the first node.
!
! Parameters
!     a  : field component in the x-direction (v or a)
!     b  : field component in the y-direction (v or a)
!     c  : field component in the z-direction (v or a)
!     d  : nonlinear term in the x-direction
!     e  : nonlinear term in the y-direction
!     f  : nonlinear term in the z-direction
!     nmb: the extension used when writting the file
!     kin: =0 computes the magnetic helicity transfer
!          =1 computes the kinetic helicity transfer
!
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1)          :: Hk,Hktot
      COMPLEX, INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX, INTENT(IN), DIMENSION(n,n,ista:iend) :: d,e,f
      COMPLEX, DIMENSION(n,n,ista:iend)             :: c1,c2,c3
      REAL                :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Hk to zero
!
      DO i = 1,n/2+1
         Hk(i) = 0.
      END DO
!
! Computes the helicity transfer
!
      tmp = 1./float(n)**6
      CALL rotor3(b,c,c1,1)
      CALL rotor3(a,c,c2,2)
      CALL rotor3(a,b,c3,3)
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn)
         DO j = 1,n
            DO k = 1,n
               kmn = int(abs(ka(k))+1)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!$omp atomic
                  Hk(kmn) = Hk(kmn)+                             &
                    (real(c1(k,j,1)*conjg(d(k,j,1)))+            &
                    real(c2(k,j,1)*conjg(e(k,j,1)))+             &
                    real(c3(k,j,1)*conjg(f(k,j,1))))*tmp
               ENDIF
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn)
            DO j = 1,n
               DO k = 1,n
                  kmn = int(abs(ka(k))+1)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!$omp atomic
                     Hk(kmn) = Hk(kmn)+                          &
                       2*(real(c1(k,j,i)*conjg(d(k,j,i)))+       &
                       real(c2(k,j,i)*conjg(e(k,j,i)))+          &
                       real(c3(k,j,i)*conjg(f(k,j,i))))*tmp
                  ENDIF
              END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn)
            DO j = 1,n
               DO k = 1,n
                  kmn = int(abs(ka(k))+1)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!$omp atomic
                     Hk(kmn) = Hk(kmn)+                          &
                       2*(real(c1(k,j,i)*conjg(d(k,j,i)))+       &
                       real(c2(k,j,i)*conjg(e(k,j,i)))+          &
                       real(c3(k,j,i)*conjg(f(k,j,i))))*tmp
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Hk,Hktot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF (kin.eq.0) THEN
            OPEN(1,file='hmtranpara.' // nmb // '.txt')
         ELSE
            OPEN(1,file='hktranpara.' // nmb // '.txt')
         ENDIF
         WRITE(1,40) Hktot
   40    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE heltpara

!*****************************************************************
      SUBROUTINE heltperp(a,b,c,d,e,f,nmb,kin)
!-----------------------------------------------------------------
!
! Computes the helicity transfer in the direction perpendicular
! to the preferred direction (rotation or uniform magnetic 
! field) in 3D Fourier space. The k-shells are cylindrical 
! surfaces (kperp=1,...,N/2+1). The output is written to a file 
! by the first node.
!
! Parameters
!     a  : field component in the x-direction (v or a)
!     b  : field component in the y-direction (v or a)
!     c  : field component in the z-direction (v or a)
!     d  : nonlinear term in the x-direction
!     e  : nonlinear term in the y-direction
!     f  : nonlinear term in the z-direction
!     nmb: the extension used when writting the file
!     kin: =0 computes the magnetic helicity transfer
!          =1 computes the kinetic helicity transfer
!
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1)            :: Hk,Hktot
      COMPLEX, INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX, INTENT(IN), DIMENSION(n,n,ista:iend) :: d,e,f
      COMPLEX, DIMENSION(n,n,ista:iend)             :: c1,c2,c3
      REAL                :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Hk to zero
!
      DO i = 1,n/2+1
         Hk(i) = 0.
      END DO
!
! Computes the helicity transfer
!
      tmp = 1./float(n)**6
      CALL rotor3(b,c,c1,1)
      CALL rotor3(a,c,c2,2)
      CALL rotor3(a,b,c3,3)
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn)
         DO j = 1,n
            kmn = int(sqrt(ka(1)**2+ka(j)**2)+.501)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
               DO k = 1,n
!$omp atomic
                  Hk(kmn) = Hk(kmn)+                             &
                    (real(c1(k,j,1)*conjg(d(k,j,1)))+            &
                    real(c2(k,j,1)*conjg(e(k,j,1)))+             &
                    real(c3(k,j,1)*conjg(f(k,j,1))))*tmp
               END DO
            ENDIF
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn)
            DO j = 1,n
               kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  DO k = 1,n
!$omp atomic
                     Hk(kmn) = Hk(kmn)+                          &
                       2*(real(c1(k,j,i)*conjg(d(k,j,i)))+       &
                       real(c2(k,j,i)*conjg(e(k,j,i)))+          &
                       real(c3(k,j,i)*conjg(f(k,j,i))))*tmp
                  END DO
               ENDIF
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn)
            DO j = 1,n
               kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  DO k = 1,n
!$omp atomic
                     Hk(kmn) = Hk(kmn)+                          &
                       2*(real(c1(k,j,i)*conjg(d(k,j,i)))+       &
                       real(c2(k,j,i)*conjg(e(k,j,i)))+          &
                       real(c3(k,j,i)*conjg(f(k,j,i))))*tmp
                  END DO
               ENDIF
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Hk,Hktot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF (kin.eq.0) THEN
            OPEN(1,file='hmtranperp.' // nmb // '.txt')
         ELSE
            OPEN(1,file='hktranperp.' // nmb // '.txt')
         ENDIF
         WRITE(1,40) Hktot
   40    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE heltperp
