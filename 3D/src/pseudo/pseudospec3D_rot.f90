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
!
! 7 Aug 2010: New specperp and 2D spectrum (T. Teitelbaum)
!=================================================================

!*****************************************************************
      SUBROUTINE specpara(a,b,c,nmb,kin,hel)
!-----------------------------------------------------------------
!
! Computes the reduced energy and helicity power spectrum 
! in the direction parallel to the preferred direction 
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
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1) :: Ek,Ektot
      DOUBLE PRECISION    :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend)             :: c1,c2,c3
      REAL(KIND=GP)       :: tmp
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
      tmp = 1.0_GP/real(n,kind=GP)**6
      IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,n
               DO k = 1,n
                  kmn = int(abs(ka(k))+1)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = (abs(a(k,j,1))**2+abs(b(k,j,1))**2+          &
                            abs(c(k,j,1))**2)*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq                       
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        tmq = 2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+     &
                                 abs(c(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        tmq = 2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+     &
                                 abs(c(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
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
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,n
               DO k = 1,n
                  kmn = int(abs(ka(k))+1)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = (abs(c1(k,j,1))**2+abs(c2(k,j,1))**2+        &
                            abs(c3(k,j,1))**2)*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        tmq = 2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+   &
                                 abs(c3(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq                          
                     ENDIF
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        tmq = 2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+   &
                              abs(c3(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
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
            Ek(i) = 0.0D0
         END DO
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,n
               DO k = 1,n
                  kmn = int(abs(ka(k))+1)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = (real(a(k,j,1)*conjg(c1(k,j,1)))+            &
                            real(b(k,j,1)*conjg(c2(k,j,1)))+            &
                            real(c(k,j,1)*conjg(c3(k,j,1))))*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+       &
                                 real(b(k,j,i)*conjg(c2(k,j,i)))+       &
                                 real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                 END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+       &
                                 real(b(k,j,i)*conjg(c2(k,j,i)))+       &
                                 real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
         CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,   &
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
! Computes the reduced energy and helicity power spectrum 
! in the direction perpendicular to the preferred direction 
! (rotation or uniform magnetic field). The k-shells are 
! cylindrical surfaces (kperp=1,...,n/2+1). It also computes 
! the spectrum of 2D modes with kz=0 for (x,y)-field 
! components and the z-field component separately. The output 
! is written to a file with three columns by the first node.
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
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1) :: Ek,Ektot
      DOUBLE PRECISION, DIMENSION(n/2+1) :: Ekp,Eptot
      DOUBLE PRECISION, DIMENSION(n/2+1) :: Ekz,Eztot
      DOUBLE PRECISION    :: tmq,tmr
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend)             :: c1,c2,c3
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: kin,hel
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.0D0
         Ekp(i) = 0.0D0
         Ekz(i) = 0.0D0
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
      tmp = 1.0_GP/real(n,kind=GP)**6
      IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq,tmr)
            DO j = 1,n
               kmn = int(sqrt(ka(1)**2+ka(j)**2)+.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  tmq = (abs(a(1,j,1))**2+abs(b(1,j,1))**2)*tmp
                  tmr = (abs(c(1,j,1))**2)*tmp
!$omp critical
                  Ekp(kmn) = Ekp(kmn)+tmq
                  Ekz(kmn) = Ekz(kmn)+tmr
                  Ek(kmn) = Ek(kmn)+tmq+tmr
!$omp end critical
                  DO k = 2,n
                     tmq = (abs(a(k,j,1))**2+abs(b(k,j,1))**2+          &
                            abs(c(k,j,1))**2)*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq,tmr)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq,tmr)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = 2*(abs(a(1,j,i))**2+abs(b(1,j,i))**2)*tmp
                     tmr = 2*abs(c(1,j,i))**2*tmp
!$omp critical
                     Ekp(kmn) = Ekp(kmn)+tmq
                     Ekz(kmn) = Ekz(kmn)+tmr
                     Ek(kmn) = Ek(kmn)+tmq+tmr
!$omp end critical
                     DO k = 2,n
                        tmq = 2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+     &
                                 abs(c(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  ENDIF
               END DO
            END DO
          ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq,tmr)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq,tmr)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = 2*(abs(a(1,j,i))**2+abs(b(1,j,i))**2)*tmp
                     tmr = 2*abs(c(1,j,i))**2*tmp
!$omp critical
                     Ekp(kmn) = Ekp(kmn)+tmq
                     Ekz(kmn) = Ekz(kmn)+tmr
                     Ek(kmn) = Ek(kmn)+tmq+tmr
!$omp end critical
                     DO k = 2,n
                        tmq = 2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+     &
                                 abs(c(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
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
!$omp parallel do private (k,kmn,tmq,tmr)
            DO j = 1,n
               kmn = int(sqrt(ka(1)**2+ka(j)**2)+.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  tmq = (abs(c1(1,j,1))**2+abs(c2(1,j,1))**2)*tmp
                  tmr = abs(c3(1,j,1))**2*tmp
!$omp critical
                  Ekp(kmn) = Ekp(kmn)+tmq
                  Ekz(kmn) = Ekz(kmn)+tmr
                  Ek(kmn) = Ek(kmn)+tmq+tmr
!$omp end critical
                  DO k = 2,n
                     tmq = (abs(c1(k,j,1))**2+abs(c2(k,j,1))**2+        &
                            abs(c3(k,j,1))**2)*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq,tmr)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq,tmr)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = 2*(abs(c1(1,j,i))**2+abs(c2(1,j,i))**2)*tmp
                     tmr = 2*abs(c3(1,j,i))**2*tmp
!$omp critical
                     Ekp(kmn) = Ekp(kmn)+tmq
                     Ekz(kmn) = Ekz(kmn)+tmr
                     Ek(kmn) = Ek(kmn)+tmq+tmr
!$omp end critical
                     DO k = 2,n
                        tmq = 2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+   &
                                 abs(c3(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq,tmr)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq,tmr)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = 2*(abs(c1(1,j,i))**2+abs(c2(1,j,i))**2)*tmp
                     tmr = 2*abs(c3(1,j,i))**2*tmp
!$omp critical
                     Ekp(kmn) = Ekp(kmn)+tmq
                     Ekz(kmn) = Ekz(kmn)+tmr
                     Ek(kmn) = Ek(kmn)+tmq+tmr
!$omp end critical
                     DO k = 2,n
                        tmq = 2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+   &
                                 abs(c3(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
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
         CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0,  &
                         MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(Ekp,Eptot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(Ekz,Eztot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) THEN
            IF (kin.eq.1) THEN
               OPEN(1,file='kspecperp.' // nmb // '.txt')
            ELSE
               OPEN(1,file='mspecperp.' // nmb // '.txt')
            ENDIF
            DO j =1,n/2+1
               WRITE(1,FMT='(E23.15,E23.15,E23.15)') Ektot(j),           &
                     Eptot(j),Eztot(j)
            END DO
            CLOSE(1)
         ENDIF
      END IF
!
! Computes the helicity spectrum
!
      IF (hel.eq.1) THEN
         DO i = 1,n/2+1
            Ek(i) = 0.0D0
            Ekp(i) = 0.0D0
            Ekz(i) = 0.0D0
         END DO
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq,tmr)
            DO j = 1,n
               kmn = int(sqrt(ka(1)**2+ka(j)**2)+.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  tmq = (real(a(1,j,1)*conjg(c1(1,j,1)))+              &
                         real(b(1,j,1)*conjg(c2(1,j,1))))*tmp
                  tmr = real(c(1,j,1)*conjg(c3(1,j,1)))*tmp
!$omp critical
                  Ekp(kmn) = Ekp(kmn)+tmq
                  Ekz(kmn) = Ekz(kmn)+tmr
                  Ek(kmn) = Ek(kmn)+tmq+tmr
!$omp end critical
                  DO k = 2,n
                     tmq = (real(a(k,j,1)*conjg(c1(k,j,1)))+           &
                            real(b(k,j,1)*conjg(c2(k,j,1)))+           &
                            real(c(k,j,1)*conjg(c3(k,j,1))))*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq,tmr)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq,tmr)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = 2*(real(a(1,j,i)*conjg(c1(1,j,i)))+         &
                           real(b(1,j,i)*conjg(c2(1,j,i))))*tmp
                     tmr = 2*real(c(1,j,i)*conjg(c3(1,j,i)))*tmp
!$omp critical
                     Ekp(kmn) = Ekp(kmn)+tmq
                     Ekz(kmn) = Ekz(kmn)+tmr
                     Ek(kmn) = Ek(kmn)+tmq+tmr
!$omp end critical
                     DO k = 2,n
                        tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+      &
                                 real(b(k,j,i)*conjg(c2(k,j,i)))+      &
                                 real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq,tmr)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq,tmr)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = 2*(real(a(1,j,i)*conjg(c1(1,j,i)))+         &
                           real(b(1,j,i)*conjg(c2(1,j,i))))*tmp
                     tmr = 2*real(c(1,j,i)*conjg(c3(1,j,i)))*tmp
!$omp critical
                     Ekp(kmn) = Ekp(kmn)+tmq
                     Ekz(kmn) = Ekz(kmn)+tmr
                     Ek(kmn) = Ek(kmn)+tmq+tmr
!$omp end critical
                     DO k = 2,n
                        tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+      &
                                 real(b(k,j,i)*conjg(c2(k,j,i)))+      &
                                 real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  ENDIF
               END DO
            END DO
         ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
         CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,  &
                         0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(Ekp,Eptot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(Ekz,Eztot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         0,MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) THEN
            IF (kin.eq.1) THEN
               OPEN(1,file='kheliperp.' // nmb // '.txt')
            ELSE IF (kin.eq.0) THEN
               OPEN(1,file='mheliperp.' // nmb // '.txt')
            ELSE
               OPEN(1,file='gheliperp.' // nmb // '.txt')
            ENDIF
            DO j =1,n/2+1
               WRITE(1,FMT='(E23.15,E23.15,E23.15)') Ektot(j),         &
                     Eptot(j),Eztot(j)
            END DO
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE specperp

!****************************************************************
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
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1) :: Ek,Ektot
      DOUBLE PRECISION    :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: d,e,f
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.0D0
      END DO
!
! Computes the kinetic energy transfer
!
      tmp = 1.0_GP/real(n,kind=GP)**6
      IF ((kin.eq.1).or.(kin.eq.2)) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,n
               DO k = 1,n
                  kmn = int(abs(ka(k))+1)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = (real(a(k,j,1)*conjg(d(k,j,1)))+         &
                            real(b(k,j,1)*conjg(e(k,j,1)))+         &
                            real(c(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(d(k,j,i)))+    &
                                 real(b(k,j,i)*conjg(e(k,j,i)))+    &
                                 real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                 END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(d(k,j,i)))+    &
                                 real(b(k,j,i)*conjg(e(k,j,i)))+    &
                                 real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
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
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,n
               DO k = 1,n
                  kmn = int(abs(ka(k))+1)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = ka2(k,j,1) *                             &
                           (real(a(k,j,1)*conjg(d(k,j,1)))+         &
                            real(b(k,j,1)*conjg(e(k,j,1)))+         &
                            real(c(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        tmq = 2*ka2(k,j,i) *                        &
                              (real(a(k,j,i)*conjg(d(k,j,i)))+      &
                               real(b(k,j,i)*conjg(e(k,j,i)))+      &
                               real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                 END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,n
                  DO k = 1,n
                     kmn = int(abs(ka(k))+1)
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        tmq = 2*ka2(k,j,i)*                         &
                              (real(a(k,j,i)*conjg(d(k,j,i)))+      &
                               real(b(k,j,i)*conjg(e(k,j,i)))+      &
                               real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
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
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1) :: Ek,Ektot
      DOUBLE PRECISION    :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: d,e,f
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.0D0
      END DO
!
! Computes the kinetic energy transfer
!
      tmp = 1.0_GP/real(n,kind=GP)**6
      IF ((kin.eq.1).or.(kin.eq.2)) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,n
               kmn = int(sqrt(ka(1)**2+ka(j)**2)+.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  DO k = 1,n
                     tmq = (real(a(k,j,1)*conjg(d(k,j,1)))+         &
                            real(b(k,j,1)*conjg(e(k,j,1)))+         &
                            real(c(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
                        tmq = 2*(real(a(k,j,i)*conjg(d(k,j,i)))+    &
                                 real(b(k,j,i)*conjg(e(k,j,i)))+    &
                                 real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
                        tmq = 2*(real(a(k,j,i)*conjg(d(k,j,i)))+    &
                                 real(b(k,j,i)*conjg(e(k,j,i)))+    &
                                 real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
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
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,n
               kmn = int(sqrt(ka(1)**2+ka(j)**2)+.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  DO k = 1,n
                     tmq = ka2(k,j,1) *                             &
                           (real(a(k,j,1)*conjg(d(k,j,1)))+         &
                            real(b(k,j,1)*conjg(e(k,j,1)))+         &
                            real(c(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
                        tmq = 2*ka2(k,j,i) *                        &
                              (real(a(k,j,i)*conjg(d(k,j,i)))+      &
                               real(b(k,j,i)*conjg(e(k,j,i)))+      &
                               real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
                        tmq = 2*ka2(k,j,i)*                         &
                              (real(a(k,j,i)*conjg(d(k,j,i)))+      &
                               real(b(k,j,i)*conjg(e(k,j,i)))+      &
                               real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
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
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1) :: Hk,Hktot
      DOUBLE PRECISION    :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend)             :: c1,c2,c3
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Hk to zero
!
      DO i = 1,n/2+1
         Hk(i) = 0.0D0
      END DO
!
! Computes the helicity transfer
!
      tmp = 1.0_GP/real(n,kind=GP)**6
      CALL rotor3(b,c,c1,1)
      CALL rotor3(a,c,c2,2)
      CALL rotor3(a,b,c3,3)
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,n
            DO k = 1,n
               kmn = int(abs(ka(k))+1)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  tmq = (real(c1(k,j,1)*conjg(d(k,j,1)))+          &
                         real(c2(k,j,1)*conjg(e(k,j,1)))+          &
                         real(c3(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                  Hk(kmn) = Hk(kmn)+tmq
               ENDIF
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
            DO j = 1,n
               DO k = 1,n
                  kmn = int(abs(ka(k))+1)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = 2*(real(c1(k,j,i)*conjg(d(k,j,i)))+     &
                              real(c2(k,j,i)*conjg(e(k,j,i)))+     &
                              real(c3(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                     Hk(kmn) = Hk(kmn)+tmq
                  ENDIF
              END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
            DO j = 1,n
               DO k = 1,n
                  kmn = int(abs(ka(k))+1)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = 2*(real(c1(k,j,i)*conjg(d(k,j,i)))+     &
                              real(c2(k,j,i)*conjg(e(k,j,i)))+     &
                              real(c3(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                     Hk(kmn) = Hk(kmn)+tmq
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
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1) :: Hk,Hktot
      DOUBLE PRECISION    :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend)             :: c1,c2,c3
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Hk to zero
!
      DO i = 1,n/2+1
         Hk(i) = 0.0D0
      END DO
!
! Computes the helicity transfer
!
      tmp = 1.0_GP/real(n,kind=GP)**6
      CALL rotor3(b,c,c1,1)
      CALL rotor3(a,c,c2,2)
      CALL rotor3(a,b,c3,3)
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,n
            kmn = int(sqrt(ka(1)**2+ka(j)**2)+.501)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
               DO k = 1,n
                  tmq = (real(c1(k,j,1)*conjg(d(k,j,1)))+          &
                         real(c2(k,j,1)*conjg(e(k,j,1)))+          &
                         real(c3(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                  Hk(kmn) = Hk(kmn)+tmq
               END DO
            ENDIF
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
            DO j = 1,n
               kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  DO k = 1,n
                     tmq = 2*(real(c1(k,j,i)*conjg(d(k,j,i)))+     &
                              real(c2(k,j,i)*conjg(e(k,j,i)))+     &
                              real(c3(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                     Hk(kmn) = Hk(kmn)+tmq
                  END DO
               ENDIF
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
            DO j = 1,n
               kmn = int(sqrt(ka(i)**2+ka(j)**2)+.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  DO k = 1,n
                     tmq = 2*(real(c1(k,j,i)*conjg(d(k,j,i)))+     &
                              real(c2(k,j,i)*conjg(e(k,j,i)))+     &
                              real(c3(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                     Hk(kmn) = Hk(kmn)+tmq                       
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

!*****************************************************************
      SUBROUTINE spec2D(a,b,c,nmb,dir,kin,hel)
!-----------------------------------------------------------------
!
! Computes the axysimmetric energy and helicity power spectrum. 
! The spectrum is angle-averaged in the azimuthal direction, 
! and depends on two wavenumbers, kperp=0,...,n/2 and 
! kpara=0,....,n/2. The output is written to a binary file by 
! the first node.
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     nmb: the extension used when writting the file
!     dir: directory where the files are written
!     kin: =2 skips energy spectrum computation
!          =1 computes the kinetic spectrum
!          =0 computes the magnetic spectrum
!     hel: =1 computes the helicity spectrum
!          =0 skips helicity spectrum computation
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend)             :: c1,c2,c3
      REAL(KIND=GP), DIMENSION(n/2+1,n/2+1)                  :: Ek,Ektot
      REAL(KIND=GP)       :: tmq,tmp
      INTEGER, INTENT(IN) :: kin,hel
      INTEGER             :: i,j,k
      INTEGER             :: kmn,kz
      CHARACTER(len=100), INTENT(IN) :: dir
      CHARACTER(len=*), INTENT(IN)   :: nmb

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         DO j = 1,n/2+1
            Ek(i,j) = 0.0_GP
         END DO
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
      tmp = 1.0_GP/real(n,kind=GP)**6
      IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kz,kmn,tmq)
            DO j = 1,n
               kmn = int(sqrt(ka(1)**2+ka(j)**2)+1.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  DO k = 1,n
                     kz = int(abs(ka(k))+1)
                     IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN
                     tmq = (abs(a(k,j,1))**2+abs(b(k,j,1))**2+        &
                            abs(c(k,j,1))**2)*tmp
!$omp atomic
                     Ek(kmn,kz) = Ek(kmn,kz)+tmq
                     ENDIF
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kz,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kz,kmn,tmq)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
                        kz = int(abs(ka(k))+1)
                        IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN
                        tmq = 2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+   &
                                 abs(c(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn,kz) = Ek(kmn,kz)+tmq
                        ENDIF
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kz,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kz,kmn,tmq)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
                        kz = int(abs(ka(k))+1)
                        IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN
                        tmq = 2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+   &
                                 abs(c(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn,kz) = Ek(kmn,kz)+tmq
                        ENDIF
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
!$omp parallel do private (k,kz,kmn,tmq)
            DO j = 1,n
               kmn = int(sqrt(ka(1)**2+ka(j)**2)+1.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  DO k = 1,n
                     kz = int(abs(ka(k))+1)
                     IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN
                     tmq = (abs(c1(k,j,1))**2+abs(c2(k,j,1))**2+      &
                            abs(c3(k,j,1))**2)*tmp
!$omp atomic
                     Ek(kmn,kz) = Ek(kmn,kz)+tmq
                     ENDIF
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kz,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kz,kmn,tmq)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
                        kz = int(abs(ka(k))+1)
                        IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN
                        tmq = 2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+ &
                                 abs(c3(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn,kz) = Ek(kmn,kz)+tmq
                        ENDIF
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kz,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kz,kmn,tmq)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
                        kz = int(abs(ka(k))+1)
                        IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN
                        tmq = 2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+ &
                                 abs(c3(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn,kz) = Ek(kmn,kz)+tmq
                        ENDIF
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
         CALL MPI_REDUCE(Ek,Ektot,(n/2+1)*(n/2+1),GC_REAL,            &
                         MPI_SUM,0,MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) THEN
            IF (kin.eq.1) THEN
               OPEN(1,file=trim(dir) // '/' // 'kspec2D.' // nmb //   &
                    '.out',form='unformatted')
            ELSE
               OPEN(1,file=trim(dir) // '/' // 'mspec2D.' // nmb //   &
                    '.out',form='unformatted')
            ENDIF
            WRITE(1) Ektot
            CLOSE(1)
         ENDIF
      ENDIF
!
! Computes the helicity spectrum
!
      IF (hel.eq.1) THEN
         DO i = 1,n/2+1
            DO j = 1,n/2+1
               Ek(i,j) = 0.0_GP
            END DO
         END DO
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kz,kmn,tmq)
            DO j = 1,n
               kmn = int(sqrt(ka(1)**2+ka(j)**2)+1.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  DO k = 1,n
                  kz = int(abs(ka(k))+1)
                  IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN
                  tmq = (real(a(k,j,1)*conjg(c1(1,j,1)))+       &
                         real(b(k,j,1)*conjg(c2(1,j,1)))+       &
                         real(c(k,j,1)*conjg(c3(k,j,1))))*tmp
!$omp atomic
                  Ek(kmn,kz) = Ek(kmn,kz)+tmq
                  ENDIF
                  ENDDO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kz,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kz,kmn,tmq)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
                        kz = int(abs(ka(k))+1)
                        IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+    &
                                 real(b(k,j,i)*conjg(c2(k,j,i)))+    &
                                 real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn,kz) = Ek(kmn,kz)+tmq
                        ENDIF
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kz,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kz,kmn,tmq)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
                        kz = int(abs(ka(k))+1)
                        IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+    &
                                 real(b(k,j,i)*conjg(c2(k,j,i)))+    &
                                 real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn,kz) = Ek(kmn,kz)+tmq
                        ENDIF
                     END DO
                  ENDIF
               END DO
            END DO
         ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
         CALL MPI_REDUCE(Ek,Ektot,(n/2+1)*(n/2+1),GC_REAL,           &
                         MPI_SUM,0,MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) THEN
            IF (kin.eq.1) THEN
               OPEN(1,file=trim(dir) // '/' // 'kheli2D.' // nmb //  &
                    '.out',form='unformatted')
            ELSE IF (kin.eq.0) THEN
               OPEN(1,file=trim(dir) // '/' // 'mheli2D.' // nmb //  &
                    '.out',form='unformatted')
            ELSE
               OPEN(1,file=trim(dir) // '/' // 'gheli2D.' // nmb //  &
                    '.out',form='unformatted')
            ENDIF
            WRITE(1) Ektot
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE spec2d
