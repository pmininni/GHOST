!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Extra subroutines to compute the passive scalar spectrum, 
! transfer function, and associated global quantities in the 
! HD, MHD, and Hall-MHD equations when a passive scalar is 
! present. You should use the FFTPLANS and MPIVARS modules 
! (see the file 'fftp_mod.f90') in each program that calls 
! any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2011 D. Rosenberg
!      NCAR
!=================================================================

!*****************************************************************
      SUBROUTINE variance(a,b,kin)
!-----------------------------------------------------------------
!
! Computes the mean variance of the passive scalar.
! The output is only valid in the first node.
!
! Parameters
!     a  : input matrix with the scalar
!     b  : at the output contains the variance
!     kin: =0 computes the variance of k^2 times the scalar
!          =1 computes the variance of the scalar
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a
      DOUBLE PRECISION, INTENT(OUT) :: b
      DOUBLE PRECISION              :: bloc
      REAL(KIND=GP)                 :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j

      bloc = 0.
      tmp = 1.0_GP/real(n,kind=GP)**4

!
! Computes the variance
!
      IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
            DO j = 1,n
                bloc = bloc+tmp*abs(a(j,1))**2
            END DO
            DO i = 2,iend
               DO j = 1,n
                   bloc = bloc+2*tmp*abs(a(j,i))**2
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                   bloc = bloc+2*tmp*abs(a(j,i))**2
               END DO
            END DO
         ENDIF
!
! Computes the variance of k^2 times the scalar
!
      ELSE IF (kin.eq.0) THEN
         IF (ista.eq.1) THEN
            DO j = 1,n
                bloc = bloc+tmp*ka2(j,1)*abs(a(j,1))**2
            END DO
            DO i = 2,iend
               DO j = 1,n
                   bloc = bloc+2*tmp*ka2(j,i)*abs(a(j,i))**2
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                   bloc = bloc+2*tmp*ka2(j,i)*abs(a(j,i))**2
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(bloc,b,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE variance

!*****************************************************************
      SUBROUTINE product(a,b,c)
!-----------------------------------------------------------------
!
! Computes the integral of the product of two scalars. 
! The output is only valid in the first node.
!
! Parameters
!     a  : first scalar
!     b  : second scalar
!     c  : at the output contains the product
!
      USE fprecision
      USE commtypes
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a,b
      DOUBLE PRECISION, INTENT(OUT) :: c
      DOUBLE PRECISION              :: cloc
      REAL(KIND=GP)                :: tmp
      INTEGER                      :: i,j,k

      cloc = 0.
      tmp = 1.0_GP/real(n,kind=GP)**4
!
! Computes the averaged inner product between the fields
!
      IF (ista.eq.1) THEN
         DO j = 1,n
             cloc = cloc+tmp*real(a(j,1)*conjg(b(j,1)))
         END DO
         DO i = 2,iend
            DO j = 1,n
                cloc = cloc+2*tmp*real(a(j,i)*conjg(b(j,i)))
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
                cloc = cloc+2*tmp*real(a(j,i)*conjg(b(j,i)))
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(cloc,c,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE product

!*****************************************************************
      SUBROUTINE pscheck(a,b,t)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of energy, 
! helicity, and null divergency of the velocity field
!
! Parameters
!     a : passive scalar concentration
!     b : source of the passive scalar
!     t : number of time steps made
!     dt: time step
!
      USE fprecision
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a,b
      DOUBLE PRECISION    :: eng,ens,pot
      REAL(KIND=GP), INTENT(IN)    :: t
      INTEGER                      :: i,j

!
! Computes the variance and the variance of k^2 times the scalar
!
      CALL variance(a,eng,1)
      CALL variance(a,ens,0)
!
! Computes the scalar injection rate
!
      CALL product(a,b,pot)
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='scalar.txt',position='append')
         WRITE(1,10) t,eng,ens,pot
   10    FORMAT( E13.6,E22.14,E22.14,E22.14 )
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE pscheck

!*****************************************************************
      SUBROUTINE spectrsc(a,nmb)
!-----------------------------------------------------------------
!
! Computes the passive scalar power spectrum. The 
! output is written to a file by the first node.
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
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1)            :: Ek,Ektot
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j
      INTEGER :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.
      END DO
!
! Computes the power spectrum
!
      tmp = 1.0_GP/real(n,kind=GP)**4
      IF (ista.eq.1) THEN
         DO j = 1,n
             kmn = int(sqrt(ka2(j,1))+0.5)
             IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                Ek(kmn) = Ek(kmn)+tmp*abs(a(j,1))**2
             ENDIF
         END DO
         DO i = 2,iend
            DO j = 1,n
                kmn = int(sqrt(ka2(j,i))+0.5)
                IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                   Ek(kmn) = Ek(kmn)+2*tmp*abs(a(j,i))**2
                ENDIF
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
                kmn = int(sqrt(ka2(j,i))+0.5)
                IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                   Ek(kmn) = Ek(kmn)+2*tmp*abs(a(j,i))**2
                ENDIF
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Exports the spectrum to a file
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='sspectrum.' // nmb // '.txt')
         WRITE(1,20) Ektot
   20    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE spectrsc

!*****************************************************************
      SUBROUTINE sctrans(a,b,nmb)
!-----------------------------------------------------------------
!
! Computes the passive scalar transfer in Fourier 
! space in 3D. The output is written to a file 
! by the first node.
!
! Parameters
!     a  : passive scalar
!     b  : nonlinear term
!     nmb: the extension used when writting the file
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1) :: Ek,Ektot
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a,b
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k
      INTEGER :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.
      END DO
!
! Computes the passive scalar transfer
!
      tmp = 1.0_GP/real(n,kind=GP)**4
      IF (ista.eq.1) THEN
         DO j = 1,n
             kmn = int(sqrt(ka2(j,1))+.501)
             IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                Ek(kmn) = Ek(kmn)+tmp*real(a(j,1)*conjg(b(j,1)))
             ENDIF
         END DO
         DO i = 2,iend
            DO j = 1,n
                kmn = int(sqrt(ka2(j,i))+.501)
                IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                   Ek(kmn) = Ek(kmn)+2*tmp*real(a(j,i)*conjg(b(j,i)))
                ENDIF
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
                kmn = int(sqrt(ka2(j,i))+.501)
                IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                   Ek(kmn) = Ek(kmn)+2*tmp*real(a(j,i)*conjg(b(j,i)))
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
         OPEN(1,file='stransfer.' // nmb // '.txt')
         WRITE(1,30) Ektot
   30    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE sctrans

!*****************************************************************
      SUBROUTINE maxabs2C(c1,c2,rmax,ic,jc)
!-----------------------------------------------------------------
!
! Computes the maximum absolute value of 
! sqrt(IFFT(c1)^2 + IFFT(c2)^2)
!
! Parameters
!     c1  : input complex field component 1
!     c2  : input complex field component 2
!     rmax: at the output contains the maximum value
!     ic  : at the output contains the i-index where the maximum is
!     jc  : at the output contains the j-index where the maximum is
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: c1,c2
      REAL(KIND=GP), DIMENSION(n,jsta:jend)    :: r1,r2
      REAL(KIND=GP), INTENT(OUT)               :: rmax
      REAL(KIND=GP)        :: bloc
      INTEGER, INTENT(OUT) :: ic,jc
      INTEGER              :: i,j

      ic = 0
      jc = 0
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      bloc = 0.0_GP
      DO j = jsta,jend
         DO i = 1,n
            bloc = max(bloc,r1(i,j)**2 + r2(i,j)**2)
         END DO
      END DO
      DO j = jsta,jend
         DO i = 1,n
            bloc = max(bloc,r1(i,j)**2+r2(i,j)**2)
         END DO
      END DO
      bloc = sqrt(bloc)/real(n,kind=GP)**2
      CALL MPI_REDUCE(bloc,rmax,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE maxabs2C


