!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines to compute spatial derivatives and nonlinear 
! terms in incompressible HD, MHD and Hall-MHD equations in 2D 
! using a pseudo-spectral method. You should use the FFTPLANS 
! and MPIVARS modules (see the file 'fftp2D_mod.f90') in each 
! program that call any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!
! 2003 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!=================================================================

!*****************************************************************
      SUBROUTINE derivk2(a,b,dir)
!-----------------------------------------------------------------
!
! Two-dimensional derivative of the matrix 'a'
!
! Parameters
!     a  : input matrix
!     b  : at the output contains the derivative da/dk_dir
!     dir: =1 derivative in the x-direction
!          =2 derivative in the y-direction
!
      USE fprecision
      USE mpivars
      USE grid
      USE kes
      USE var
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(n,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,ista:iend) :: b
      INTEGER, INTENT(IN) :: dir
      INTEGER :: i,j

!
! Derivative in the x-direction
!
      IF (dir.eq.1) THEN
         DO i = ista,iend
            DO j = 1,n
               b(j,i) = im*ka(i)*a(j,i)
            END DO
         END DO
!
! Derivative in the y-direction
!
      ELSE
         DO i = ista,iend
            DO j = 1,n
               b(j,i) = im*ka(j)*a(j,i)
            END DO
         END DO
      ENDIF

      RETURN
      END SUBROUTINE derivk2

!*****************************************************************
      SUBROUTINE laplak2(a,b)
!-----------------------------------------------------------------
!
! Two-dimensional Laplacian of the matrix 'a'
!
! Parameters
!     a: input matrix
!     b: at the output contains the Laplacian d2a/dka2
!
      USE fprecision
      USE mpivars
      USE grid
      USE kes
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(n,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,ista:iend) :: b
      INTEGER :: i,j

      DO i = ista,iend
         DO j = 1,n
            b(j,i) = -ka2(j,i)*a(j,i)
         END DO
      END DO

      RETURN
      END SUBROUTINE laplak2

!*****************************************************************
      SUBROUTINE poisson(a,b,c)
!-----------------------------------------------------------------
!
! Poisson bracket of the scalar fields A and B 
! in real space.
!
! Parameters
!     a: input matrix
!     b: input matrix
!     c: Poisson bracket {A,B} [output]
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(n,ista:iend) :: a,b
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,ista:iend) :: c
      COMPLEX(KIND=GP), DIMENSION(n,ista:iend) :: c1,c2
      REAL(KIND=GP), DIMENSION(n,jsta:jend)    :: r1,r2,r3
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j

!
! Computes dA/dx.dB/dy
!
      CALL derivk2(a,c1,1)
      CALL derivk2(b,c2,2)
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,n
            r3(i,j) = r1(i,j)*r2(i,j)
         END DO
      END DO
!
! Computes dA/dy.dB/dx
!
      CALL derivk2(a,c1,2)
      CALL derivk2(b,c2,1)
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      tmp = 1.0_GP/real(n,kind=GP)**4
      DO j = jsta,jend
         DO i = 1,n
            r3(i,j) = (r3(i,j)-r1(i,j)*r2(i,j))*tmp
         END DO
      END DO

      CALL fftp2d_real_to_complex(planrc,r3,c,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE poisson

!*****************************************************************
      SUBROUTINE energy(a,b,kin)
!-----------------------------------------------------------------
!
! Computes the mean kinetic or magnetic energy in 2D,
! and the mean square current density or vorticity. 
! The output is valid only in the first node.
!
! Parameters
!     a  : input matrix with the scalar field
!     b  : at the output contains the energy
!     kin: =2 computes the square of the scalar field
!          =1 computes the energy
!          =0 computes the current or vorticity
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a
      DOUBLE PRECISION, INTENT(OUT) :: b
      DOUBLE PRECISION              :: bloc
      REAL(KIND=GP)                 :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j

      bloc = 0.0D0
      tmp = 1.0_GP/real(n,kind=GP)**4
!
! Computes the square of the scalar field
!
      IF (kin.eq.2) THEN
         IF (ista.eq.1) THEN
            DO j = 1,n
               bloc = bloc+abs(a(j,1))**2*tmp
            END DO
            DO i = 2,iend
               DO j = 1,n
                  bloc = bloc+2*abs(a(j,i))**2*tmp
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  bloc = bloc+2*abs(a(j,i))**2*tmp
               END DO
            END DO
         ENDIF
!
! Computes the energy
!
      ELSE IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
            DO j = 1,n
               bloc = bloc+ka2(j,1)*abs(a(j,1))**2*tmp
            END DO
            DO i = 2,iend
               DO j = 1,n
                  bloc = bloc+2*ka2(j,i)*abs(a(j,i))**2*tmp
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  bloc = bloc+2*ka2(j,i)*abs(a(j,i))**2*tmp
               END DO
            END DO
         ENDIF
!
! Computes the current or vorticity
!
      ELSE
         IF (ista.eq.1) THEN
            DO j = 1,n
               bloc = bloc+ka2(j,1)**2*abs(a(j,1))**2*tmp
            END DO
            DO i = 2,iend
               DO j = 1,n
                  bloc = bloc+2*ka2(j,i)**2*abs(a(j,i))**2*tmp
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  bloc = bloc+2*ka2(j,i)**2*abs(a(j,i))**2*tmp
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
      END SUBROUTINE energy

!*****************************************************************
      SUBROUTINE maxabs(a,b,c,d)
!-----------------------------------------------------------------
!
! Computes the maximum absolute value of the 
! vorticity or current density. The output is 
! only valid in the first node.
!
! Parameters
!     a: input field
!     b: at the output contains the maximum value
!     c: at the output contains the i-index where the maximum is
!     d: at the output contains the j-index where the maximum is
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a
      COMPLEX(KIND=GP), DIMENSION(n,ista:iend) :: c1
      REAL(KIND=GP), DIMENSION(n,jsta:jend)    :: r1
      REAL(KIND=GP), INTENT(OUT)               :: b
      REAL(KIND=GP)        :: bloc
      INTEGER, INTENT(OUT) :: c,d
      INTEGER              :: i,j

      CALL laplak2(a,c1)
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      bloc = 0.0_GP
      DO j = jsta,jend
         DO i = 1,n
            bloc = max(bloc,abs(r1(i,j)))
         END DO
      END DO
      bloc = bloc/real(n,kind=GP)**2
      CALL MPI_REDUCE(bloc,b,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE maxabs

!*****************************************************************
      SUBROUTINE hdcheck(a,b,t)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of energy in HD 2D
!
! Parameters
!     a  : streamfunction
!     b  : external force
!     t  : time
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a,b
      DOUBLE PRECISION :: eng,ens,pot,tmp
      REAL(KIND=GP) :: t
      REAL(KIND=GP) :: tmq
      INTEGER       :: i,j

      tmq = 1.0_GP/real(n,kind=GP)**4
!
! Computes the mean energy and enstrophy
!
      CALL energy(a,eng,1)
      CALL energy(a,ens,0)
!
! Computes the energy injection rate
!
      tmp = 0.0D0
      IF (ista.eq.1) THEN
         DO j = 1,n
            tmp = tmp+ka2(j,1)*real(b(j,1)*conjg(a(j,1)))*tmq
         END DO
         DO i = 2,iend
            DO j = 1,n
               tmp = tmp+2*ka2(j,i)*real(b(j,i)*conjg(a(j,i)))*tmq
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               tmp = tmp+2*ka2(j,i)*real(b(j,i)*conjg(a(j,i)))*tmq
            END DO
         END DO
      ENDIF
      CALL MPI_REDUCE(tmp,pot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='balance.txt',position='append')
         WRITE(1,10) t,eng,ens,pot
   10    FORMAT( E13.6,E26.18,E26.18,E26.18 )
         CLOSE(1)
      ENDIF      

      RETURN
      END SUBROUTINE hdcheck

!*****************************************************************
      SUBROUTINE spectrum(a,nmb,kin)
!-----------------------------------------------------------------
!
! Computes the energy power spectrum in 2D. 
! The output is written to a file by the first node.
!
! Parameters
!     a  : streamfunction or vector potential
!     nmb: the extension used when writting the file
!     kin: =0 computes the magnetic spectrum
!          =1 computes the kinetic spectrum
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE filefmt
      USE grid
      USE kes
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a
      DOUBLE PRECISION, DIMENSION(n/2+1) :: Ek,Ektot
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: kmn
      INTEGER             :: i,j
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.0D0
      END DO
!
! Computes the energy spectrum
!
      tmp = 1.0_GP/real(n,kind=GP)**4
      IF (ista.eq.1) THEN
         DO j = 1,n
            kmn = int(sqrt(ka2(j,1))+.5)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
               Ek(kmn) = Ek(kmn)+ka2(j,1)*abs(a(j,1))**2*tmp
            ENDIF
         END DO
         DO i = 2,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  Ek(kmn) = Ek(kmn)+2*ka2(j,i)*abs(a(j,i))**2*tmp
               ENDIF
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  Ek(kmn) = Ek(kmn)+2*ka2(j,i)*abs(a(j,i))**2*tmp
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
         IF (kin.eq.0) THEN
            OPEN(1,file='mspectrum.' // nmb // '.txt')
         ELSE IF (kin.eq.1) THEN
            OPEN(1,file='kspectrum.' // nmb // '.txt')
         ENDIF
         WRITE(1,20) Ektot
   20    FORMAT( E23.15 )
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE spectrum
