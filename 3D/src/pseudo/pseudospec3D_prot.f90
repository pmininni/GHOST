!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Extra subroutines to compute the variance spectrum and 
! transfer functions for a passive scalar in 3D in the 
! anisotropic case (e.g. in the rotating frame or with an 
! imposed external magnetic field). Quantities parallel 
! and perpendicular to the z direction in Fourier space 
! can be computed. You should use the FFTPLANS and MPIVARS 
! modules (see the file 'fftp_mod.f90') in each program 
! that calls any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2009 Paola Rodriguez Imazio and Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!
! 7 Aug 2010: New specperp and 2D spectrum (Teitelbaum & Mininni)
!=================================================================

!*****************************************************************
      SUBROUTINE specscpa(a,nmb)
!-----------------------------------------------------------------
!
! Computes the reduced power spectrum of the passive scalar 
! in the direction parallel to the preferred direction 
! (rotation or uniform magnetic field). As a result, the 
! k-shells are planes with normal (0,0,kz) (kz=0,...,n/2). 
! The output is written to a file by the first node.
!
! Parameters
!     a  : input matrix with the passive scalar
!     nmb: the extension used when writting the file

      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1) :: Ek,Ektot
      DOUBLE PRECISION :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a
      REAL(KIND=GP)    :: tmp
      INTEGER          :: i,j,k
      INTEGER          :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.0D0
      END DO
!
! Computes the power spectrum
!
      tmp = 1.0_GP/real(n,kind=GP)**6
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,n
            DO k = 1,n
               kmn = int(abs(ka(k))+1)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  tmq = (abs(a(k,j,1))**2)*tmp
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
                     tmq = 2*(abs(a(k,j,i))**2)*tmp
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
                     tmq = 2*(abs(a(k,j,i))**2)*tmp
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
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         OPEN(1,file='sspecpara.' // nmb // '.txt')
         WRITE(1,10) Ektot
   10    FORMAT( E23.15 )
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE specscpa

!*****************************************************************
 SUBROUTINE specscpe(a,nmb)
!-----------------------------------------------------------------
!
! Computes the reduced power spectrum of the passive scalar 
! in the direction perpendicular to the preferred direction 
! (rotation or uniform magnetic field). The k-shells are 
! cylindrical surfaces (kperp=1,...,n/2+1). It also computes 
! the spectrum of 2D modes with kz=0. The output is written 
! to a file with two columns by the first node.
!
! Parameters
!     a  : input matrix with the passive scalar
!     nmb: the extension used when writting the file
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
      DOUBLE PRECISION :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a
      REAL(KIND=GP)    :: tmp
      INTEGER          :: i,j,k
      INTEGER          :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.0D0
         Ekp(i) = 0.0D0
      END DO
!
! Computes the power spectrum
!
      tmp = 1.0_GP/real(n,kind=GP)**6
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,n
            kmn = int(sqrt(ka(1)**2+ka(j)**2)+.501)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
               tmq = (abs(a(1,j,1))**2)*tmp
!$omp critical
               Ekp(kmn) = Ekp(kmn)+tmq
               Ek(kmn) = Ek(kmn)+tmq
!$omp end critical
               DO k = 2,n
                  tmq = (abs(a(k,j,1))**2)*tmp
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
                  tmq = 2*(abs(a(1,j,i))**2)*tmp
!$omp critical
                  Ekp(kmn) = Ekp(kmn)+tmq
                  Ek(kmn) = Ek(kmn)+tmq
!$omp end critical
                  DO k = 2,n
                     tmq = 2*(abs(a(k,j,i))**2)*tmp
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
                  tmq = 2*(abs(a(1,j,i))**2)*tmp
!$omp critical
                  Ekp(kmn) = Ekp(kmn)+tmq
                  Ek(kmn) = Ek(kmn)+tmq
!$omp end critical
                  DO k = 2,n
                     tmq = 2*(abs(a(k,j,i))**2)*tmp
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
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0,  &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(Ekp,Eptot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         OPEN(1,file='sspecperp.' // nmb // '.txt')
         DO j =1,n/2+1
         WRITE(1,FMT='(E23.15,E23.15)') Ektot(j),Eptot(j)
         END DO
         CLOSE(1)
      ENDIF
!
      RETURN
      END SUBROUTINE specscpe

!*****************************************************************
      SUBROUTINE sctpara(a,b,nmb,isc)
!-----------------------------------------------------------------
!
! Computes the transfer function for the passive scalar in 
! the direction parallel to the preferred direction (rotation 
! or uniform magnetic field) in 3D Fourier space. The 
! k-shells are planes with normal (0,0,kz) (kz=0,...,n/2). 
! The output is written to a file by the first node.
!
! Parameters
!     a  : passive scalar
!     b  : nonlinear term
!     nmb: the extension used when writting the file
!     isc: scalar id, used if isc >= 0
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
      DOUBLE PRECISION :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b
      REAL(KIND=GP)    :: tmp
      INTEGER         , INTENT(IN)                           :: isc
      INTEGER          :: i,j,k
      INTEGER          :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb
      CHARACTER(len=1)             :: si

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.0D0
      END DO
!
! Computes the passive scalar transfer
!
      tmp = 1.0_GP/real(n,kind=GP)**6
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,n
            DO k = 1,n
               kmn = int(abs(ka(k))+1)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  tmq = tmp*real(a(k,j,1)*conjg(b(k,j,1)))
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
                     tmq = 2*tmp*real(a(k,j,i)*conjg(b(k,j,i)))
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
                     tmq = 2*tmp*real(a(k,j,i)*conjg(b(k,j,i)))
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
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF ( isc.GE.0 ) THEN
           WRITE(si,'(i1.1)')isc
           OPEN(1,file='s' // trim(si) // 'tranpara.' // nmb // '.txt')
         ELSE
           OPEN(1,file='stranpara.' // nmb // '.txt')
         ENDIF
         WRITE(1,30) Ektot
   30    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE sctpara

!*****************************************************************
      SUBROUTINE sctperp(a,b,nmb,isc)
!-----------------------------------------------------------------
!
! Computes the transfer function for the passive scalar in 
! the direction perpendicular to the preferred direction 
! (rotation or uniform magnetic field) in 3D Fourier space. 
! The k-shells are cylindrical surfaces (kperp=1,...,N/2+1). 
! The output is written to a file by the first node.
!
! Parameters
!     a  : passive scalar
!     b  : nonlinear term
!     nmb: the extension used when writting the file
!     isc: scalar id, used if isc >= 0
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
      DOUBLE PRECISION :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b
      REAL(KIND=GP)    :: tmp
      INTEGER         , INTENT(IN)                            :: isc
      INTEGER          :: i,j,k
      INTEGER          :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb
      CHARACTER(len=1)             :: si

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.0D0
      END DO
!
! Computes the passive scalar transfer
!
      tmp = 1.0_GP/real(n,kind=GP)**6
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,n
            kmn = int(sqrt(ka(1)**2+ka(j)**2)+.501)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
               DO k = 1,n
                  tmq = tmp*real(a(k,j,1)*conjg(b(k,j,1)))
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
                     tmq = 2*tmp*real(a(k,j,i)*conjg(b(k,j,i)))
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
                     tmq = 2*tmp*real(a(k,j,i)*conjg(b(k,j,i)))
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
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
        IF ( isc.GE.0 ) THEN
           WRITE(si,'(i1.1)')isc
           OPEN(1,file='s' // trim(si) // 'tranpara.' // nmb // '.txt')
         ELSE
           OPEN(1,file='stranperp.' // nmb // '.txt')
         ENDIF
         WRITE(1,40) Ektot
   40    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE sctperp

!*****************************************************************
      SUBROUTINE specsc2D(a,nmb,dir)
!-----------------------------------------------------------------
!
! Computes the axysimmetric power spectrum of the passive 
! scalar. The spectrum is angle-averaged in the azimuthal 
! direction, and depends on two wavenumbers, kperp=0,...,n/2 
! and kpara=0,....,n/2. The output is written to a binary file 
! by the first node.
!
! Parameters
!     a  : input matrix with the passive scalar
!     nmb: the extension used when writting the file
!     dir: directory where the files are written
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a
      REAL(KIND=GP), DIMENSION(n/2+1,n/2+1)                  :: Ek,Ektot
      REAL(KIND=GP)       :: tmq,tmp
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
! Computes the power spectrum
!
      tmp = 1.0_GP/real(n,kind=GP)**6
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kz,kmn,tmq)
         DO j = 1,n
            kmn = int(sqrt(ka(1)**2+ka(j)**2)+1.501)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
               DO k = 1,n
                  kz = int(abs(ka(k))+1)
                  IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN
                  tmq = (abs(a(k,j,1))**2)*tmp
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
                     tmq = 2*(abs(a(k,j,i))**2)*tmp
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
                     tmq = 2*(abs(a(k,j,i))**2)*tmp
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
      CALL MPI_REDUCE(Ek,Ektot,(n/2+1)*(n/2+1),GC_REAL,         &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(dir) // '/' // 'sspec2D.' // nmb //   &
              '.out',form='unformatted')
         WRITE(1) Ektot
         CLOSE(1)
      ENDIF

      END SUBROUTINE specsc2d
