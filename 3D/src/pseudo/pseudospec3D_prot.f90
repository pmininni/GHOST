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
      SUBROUTINE specscpa(a,nmb,isc)
!-----------------------------------------------------------------
!
! Computes the reduced power spectrum of the passive scalar in
! the direction parallel to the preferred direction (rotation
! or uniform magnetic field). As a result, the k-shells 
! are planes with normal (0,0,kz), kz = Dkz*(0,...,nz/2).
! Normalization of the reduced spectrum is such that
! E = sum[E(kz).Dkz], where Dkz is the width of the Fourier shells
! in kz. The output is written to a file by the first node.
!
! Output files contain:
! 'sspecpara.XXX.txt' : kz, V(kz) (power spectrum of the scalar)
! 'sNspecpara.XXX.txt': kz, V(kz) (same for the N-th scalar)
!
! Parameters
!     a  : input matrix with the passive scalar
!     nmb: the extension used when writting the file
!     isc: index to specify which scalar the spectrum 
!          represents; modifies output file name
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nz/2+1) :: Ek,Ektot
      DOUBLE PRECISION :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a
      REAL(KIND=GP)    :: tmp
      INTEGER,          INTENT(IN)                           :: isc
      INTEGER          :: i,j,k
      INTEGER          :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb
      CHARACTER(len=1)             :: si

!
! Sets Ek to zero
!
      DO k = 1,nz/2+1
         Ek(k) = 0.0D0
      END DO
!
! Computes the power spectrum
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,ny
            DO k = 1,nz
               kmn = int(abs(kz(k))*Lz+1)
               IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                  tmq = (abs(a(k,j,1))**2)*tmp
!$omp atomic
                  Ek(kmn) = Ek(kmn)+tmq
               ENDIF
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(abs(kz(k))*Lz+1)
                  IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
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
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(abs(kz(k))*Lz+1)
                  IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
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
      CALL MPI_REDUCE(Ek,Ektot,nz/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF ( isc.gt.0 ) THEN
           WRITE(si,'(i1.1)') isc
           OPEN(1,file='s' // si // 'specpara.' // nmb // '.txt')
         ELSE
           OPEN(1,file='sspecpara.' // nmb // '.txt')
         ENDIF
         DO k = 1,nz/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkz*(k-1),Ektot(k)*Lz
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE specscpa

!*****************************************************************
 SUBROUTINE specscpe(a,nmb,isc)
!-----------------------------------------------------------------
!
! Computes the reduced power spectrum of the passive scalar 
! in the direction perpendicular to the preferred direction 
! (rotation or uniform magnetic field). The k-shells are 
! cylindrical surfaces with
! kperp = Dkk*(0,...,max{nx*Dkx/Dkk,nyDky/Dkk}/2). It also
! computes the spectrum of 2D modes with kz=0.
! Normalization of the reduced spectrum is such that
! E = sum[E(kperp).Dkk], where Dkk is the width of the Fourier
! shells. The output is written to a file with two columns by
! the first node.
!
! Output files contain [kp = Dkk*sqrt(kx**2+ky**2)]:
! 'sspecperp.XXX.txt' : kp, V(kp), v(kp,kz=0)
! 'sNspecperp.XXX.txt': kp, V(kp), v(kp,kz=0) (for the N-th scalar)
!
! Parameters   
!     a  : input matrix with the passive scalar
!     nmb: the extension used when writting the file
!     isc: index to specify which scalar the spectrum 
!          represents; modifies output file name
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nmaxperp/2+1) :: Ek,Ektot
      DOUBLE PRECISION, DIMENSION(nmaxperp/2+1) :: Ekp,Eptot
      DOUBLE PRECISION :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a
      REAL(KIND=GP)    :: tmp
      INTEGER,          INTENT(IN)                           :: isc
      INTEGER          :: i,j,k
      INTEGER          :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb
      CHARACTER(len=1)             :: si

!
! Sets Ek to zero
!
      DO i = 1,nmaxperp/2+1
         Ek (i) = 0.0D0
         Ekp(i) = 0.0D0
      END DO
!
! Computes the power spectrum
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,ny
            kmn = int(sqrt(kx(1)**2+ky(j)**2)/Dkk+1)
            IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
               tmq = (abs(a(1,j,1))**2)*tmp
!$omp critical
               Ekp(kmn) = Ekp(kmn)+tmq
               Ek(kmn) = Ek(kmn)+tmq
!$omp end critical
               DO k = 2,nz
                  tmq = (abs(a(k,j,1))**2)*tmp
!$omp atomic
                  Ek(kmn) = Ek(kmn)+tmq
               END DO
            ENDIF
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
            DO j = 1,ny
               kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
               IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                  tmq = 2*(abs(a(1,j,i))**2)*tmp
!$omp critical
                  Ekp(kmn) = Ekp(kmn)+tmq
                  Ek(kmn) = Ek(kmn)+tmq
!$omp end critical
                  DO k = 2,nz
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
            DO j = 1,ny
               kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
               IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                  tmq = 2*(abs(a(1,j,i))**2)*tmp
!$omp critical
                  Ekp(kmn) = Ekp(kmn)+tmq
                  Ek(kmn) = Ek(kmn)+tmq
!$omp end critical
                  DO k = 2,nz
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
      CALL MPI_REDUCE(Ek,Ektot,nmaxperp/2+1,MPI_DOUBLE_PRECISION,  &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(Ekp,Eptot,nmaxperp/2+1,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF ( isc.gt.0 ) THEN
           WRITE(si,'(i1.1)') isc
           OPEN(1,file='s' // si // 'specperp.' // nmb // '.txt')
         ELSE
           OPEN(1,file='sspecperp.' // nmb // '.txt')
         ENDIF
         DO j = 1,nmaxperp/2+1
            WRITE(1,FMT='(E23.15,E23.15,E23.15)') Dkk*(j-1), &
                                 Ektot(j)/Dkk, Eptot(j)/Dkk
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
! or uniform magnetic field) in 3D Fourier space. The k-shells
! are planes with normal (0,0,kz), kz = Dkz*(0,...,nz/2).
! Normalization of the transfer function is such that the
! flux is Pi = -sum[T(kz).Dkz], where Dkz is the width of the
! Fourier shells. The output is written to a file by the
! first node.
!
! Output files contain:
! 'stranpara.XXX.txt' : kz, Ts(kz) (scalar transfer function)
! 'sNtranpara.XXX.txt': kz, Ts(kz) (same for the N-th scalar)
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
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nz/2+1) :: Ek,Ektot
      DOUBLE PRECISION :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      REAL(KIND=GP)    :: tmp
      INTEGER         , INTENT(IN)                             :: isc
      INTEGER          :: i,j,k
      INTEGER          :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb
      CHARACTER(len=1)             :: si

!
! Sets Ek to zero
!
      DO k = 1,nz/2+1
         Ek(k) = 0.0D0
      END DO
!
! Computes the passive scalar transfer
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,ny
            DO k = 1,nz
               kmn = int(abs(kz(k))*Lz+1)
               IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                  tmq = tmp*real(a(k,j,1)*conjg(b(k,j,1)))
!$omp atomic
                  Ek(kmn) = Ek(kmn)+tmq
               ENDIF
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(abs(kz(k))*Lz+1)
                  IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
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
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(abs(kz(k))*Lz+1)
                  IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
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
      CALL MPI_REDUCE(Ek,Ektot,nz/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF ( isc.gt.0 ) THEN
           WRITE(si,'(i1.1)')isc
           OPEN(1,file='s' // trim(si) // 'tranpara.' // nmb // '.txt')
         ELSE
           OPEN(1,file='stranpara.' // nmb // '.txt')
         ENDIF
         DO k = 1,nz/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkz*(k-1),Ektot(k)*Lz
         END DO
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
! The k-shells are cylindrical surfaces with
! kperp = Dkk*(0,...,max{nx*Dkx/Dkk,nyDky/Dkk}/2).
! Normalization of the transfer function is such that the
! flux is Pi = -sum[T(kperp).Dkk], where Dkk is the width of
! the Fourier shells. The output is written to a file by the 
! first node.
!
! Output files contain [kp = Dkk*sqrt(kx**2+ky**2)]:
! 'stranperp.XXX.txt' : kp, Ts(kp) (scalar transfer function)
! 'sNtranperp.XXX.txt': kp, Ts(kp) (same for the N-th scalar)
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
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nmaxperp/2+1) :: Ek,Ektot
      DOUBLE PRECISION :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b
      REAL(KIND=GP)    :: tmp
      INTEGER         , INTENT(IN)                            :: isc
      INTEGER          :: i,j,k
      INTEGER          :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb
      CHARACTER(len=1)             :: si

!
! Sets Ek to zero
!
      DO i = 1,nmaxperp/2+1
         Ek(i) = 0.0D0
      END DO
!
! Computes the passive scalar transfer
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,ny
            kmn = int(sqrt(kx(1)**2+ky(j)**2)/Dkk+1)
            IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
               DO k = 1,nz
                  tmq = tmp*real(a(k,j,1)*conjg(b(k,j,1)))
!$omp atomic
                  Ek(kmn) = Ek(kmn)+tmq
               END DO
            ENDIF
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
            DO j = 1,ny
               kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
               IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                  DO k = 1,nz
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
            DO j = 1,ny
               kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
               IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                  DO k = 1,nz
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
      CALL MPI_REDUCE(Ek,Ektot,nmaxperp/2+1,MPI_DOUBLE_PRECISION,     &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
        IF ( isc.gt.0 ) THEN
           WRITE(si,'(i1.1)')isc
           OPEN(1,file='s' // trim(si) // 'tranpara.' // nmb // '.txt')
         ELSE
           OPEN(1,file='stranperp.' // nmb // '.txt')
         ENDIF
         DO j = 1,nmaxperp/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(j-1),Ektot(j)/Dkk
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE sctperp

!*****************************************************************
      SUBROUTINE specsc2D(a,nmb,dir,isc)
!-----------------------------------------------------------------
!
! Computes the axysimmetric power spectrum of the passive 
! scalar. The spectrum is angle-averaged in the azimuthal 
! direction, and depends on two wavenumbers, (kperp,kpara),
! with kperp = Dkk*(0,...,max{nx*Dkx/Dkk,nyDky/Dkk}/2) and
! kpara = Dkz*(0,...,nz/2). The actual values of these
! wavenumbers can be found in the first column of 'kspecpara'
! and 'kspecperp' files. This spectrum is not normalized (i.e.,
! it is not divided by the area of the bins in Fourier space,
! nor by sin(theta) to obtain curcles in the isotropic case).
! The output is written to a binary file by the first node.
!
! Output files contain:
! 'odir/sspec2D.XXX.out' : 2D spectrum v(kperp,kpara)
! 'odir/sNspec2D.XXX.out': Same for the N-th scalar
!
! Parameters
!     a  : input matrix with the passive scalar
!     nmb: the extension used when writting the file
!     dir: directory where the files are written
!     isc: if > 0, changes filename prefix
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a
      REAL(KIND=GP), DIMENSION(nmaxperp/2+1,nz/2+1)     :: Ek,Ektot
      REAL(KIND=GP)       :: tmq,tmp
      INTEGER             :: i,j,k
      INTEGER, INTENT(IN) :: isc
      INTEGER             :: kmn,kmz
      CHARACTER(len=100), INTENT(IN) :: dir
      CHARACTER(len=*), INTENT(IN)   :: nmb
      CHARACTER(len=1)               :: si

!
! Sets Ek to zero
!
      DO i = 1,nmaxperp/2+1
         DO k = 1,nz/2+1
            Ek(i,k) = 0.0_GP
         END DO
      END DO
!
! Computes the power spectrum
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmz,kmn,tmq)
         DO j = 1,ny
            kmn = int(sqrt(kx(1)**2+ky(j)**2)/Dkk+1)
            IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
               DO k = 1,nz
                  kmz = int(abs(kz(k))*Lz+1)
                  IF ((kmz.gt.0).and.(kmz.le.nz/2+1)) THEN
                  tmq = (abs(a(k,j,1))**2)*tmp
!$omp atomic
                  Ek(kmn,kmz) = Ek(kmn,kmz)+tmq
                  ENDIF
               END DO
            ENDIF
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmz,kmn,tmq)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmz,kmn,tmq)
            DO j = 1,ny
               kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
               IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                  DO k = 1,nz
                     kmz = int(abs(kz(k))*Lz+1)
                     IF ((kmz.gt.0).and.(kmz.le.nz/2+1)) THEN
                     tmq = 2*(abs(a(k,j,i))**2)*tmp
!$omp atomic
                     Ek(kmn,kmz) = Ek(kmn,kmz)+tmq
                     ENDIF
                  END DO
               ENDIF
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmz,kmn,tmq)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmz,kmn,tmq)
            DO j = 1,ny
               kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
               IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                  DO k = 1,nz
                     kmz = int(abs(kz(k))*Lz+1)
                     IF ((kmz.gt.0).and.(kmz.le.nz/2+1)) THEN
                     tmq = 2*(abs(a(k,j,i))**2)*tmp
!$omp atomic
                     Ek(kmn,kmz) = Ek(kmn,kmz)+tmq
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
      CALL MPI_REDUCE(Ek,Ektot,(nmaxperp/2+1)*(nz/2+1),GC_REAL,       &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF ( isc.gt.0 ) THEN
           WRITE(si,'(i1.1)') isc
           OPEN(1,file=trim(dir) // '/' // si // 'spec2D.' // nmb //  &
                '.out',form='unformatted')
         ELSE
           OPEN(1,file=trim(dir) // '/' // 'sspec2D.' // nmb //       &
                '.out',form='unformatted')
         ENDIF
         WRITE(1) Ektot
         CLOSE(1)
      ENDIF

      END SUBROUTINE specsc2d
