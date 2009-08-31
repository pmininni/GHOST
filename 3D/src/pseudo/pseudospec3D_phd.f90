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
! 2009 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!=================================================================

!*****************************************************************
      SUBROUTINE advect3(a,b,c,d,e)
!-----------------------------------------------------------------
!
! Three-dimensional inner product -A.grad(B) in 
! real space. The components of the field A are 
! given by the arrays a, b and c, B is a scalar 
! quantity given by d.
!
! Parameters
!     a: input matrix in the x-direction
!     b: input matrix in the y-direction
!     c: input matrix in the z-direction
!     d: input matrix with the passive scalar
!     e: product (A.grad)B in Fourier space [output]
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend)  :: a,b,c,d
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,n,ista:iend) :: e
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend) :: c1,c2
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)    :: r1,r2
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)    :: r3
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k

!
! Computes (A_x.dx)B
!
      c1 = a
      CALL derivk3(d,c2,1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)

      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r3(i,j,k) = r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO
!
! Computes (A_y.dy)B
!
      c1 = b
      CALL derivk3(d,c2,2)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)

      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r3(i,j,k) = r3(i,j,k)+r1(i,j,k)*r2(i,j,k)
            END DO
         END DO
      END DO
!
! Computes (A_z.dz)B
!
      c1 = c
      CALL derivk3(d,c2,3)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)

      tmp = -1./real(n,kind=GP)**6   !we need -A.grad(B)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r3(i,j,k) = (r3(i,j,k)+r1(i,j,k)*r2(i,j,k))*tmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,r3,e,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE advect3

!*****************************************************************
      SUBROUTINE variance(a,b,kin)
!-----------------------------------------------------------------
!
! Computes the mean variance of the passive scalar.
! The output is only valid in the first node.
!
! Parameters
!     a  : input matrix with the scalar
!     d  : at the output contains the variance
!     kin: =0 computes the variance of k^2 times the scalar
!          =1 computes the variance of the scalar
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a
      DOUBLE PRECISION, INTENT(OUT) :: b
      DOUBLE PRECISION              :: bloc
      REAL(KIND=GP)                :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k

      bloc = 0.
      tmp = 1./real(n,kind=GP)**6

!
! Computes the variance
!
      IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
            DO j = 1,n
               DO k = 1,n
                  bloc = bloc+tmp*abs(a(k,j,1))**2
               END DO
            END DO
            DO i = 2,iend
               DO j = 1,n
                  DO k = 1,n
                     bloc = bloc+2*tmp*abs(a(k,j,i))**2
                  END DO
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  DO k = 1,n
                     bloc = bloc+2*tmp*abs(a(k,j,i))**2
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the variance of k^2 times the scalar
!
      ELSE IF (kin.eq.0) THEN
         IF (ista.eq.1) THEN
            DO j = 1,n
               DO k = 1,n
                  bloc = bloc+tmp*ka2(k,j,1)*abs(a(k,j,1))**2
               END DO
            END DO
            DO i = 2,iend
               DO j = 1,n
                  DO k = 1,n
                     bloc = bloc+2*tmp*ka2(k,j,i)*abs(a(k,j,i))**2
                  END DO
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  DO k = 1,n
                     bloc = bloc+2*tmp*ka2(k,j,i)*abs(a(k,j,i))**2
                  END DO
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

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b
      DOUBLE PRECISION, INTENT(OUT) :: c
      DOUBLE PRECISION              :: cloc
      REAL(KIND=GP)                :: tmp
      INTEGER             :: i,j,k

      cloc = 0.
      tmp = 1./real(n,kind=GP)**6
!
! Computes the averaged inner product between the fields
!
      IF (ista.eq.1) THEN
         DO j = 1,n
            DO k = 1,n
               cloc = cloc+tmp*real(a(k,j,1)*conjg(b(k,j,1)))
            END DO
         END DO
         DO i = 2,iend
            DO j = 1,n
               DO k = 1,n
                  cloc = cloc+2*tmp*real(a(k,j,i)*conjg(b(k,j,i)))
               END DO
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               DO k = 1,n
                  cloc = cloc+2*tmp*real(a(k,j,i)*conjg(b(k,j,i)))
               END DO
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
      SUBROUTINE pscheck(a,b,t,dt)
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

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b
      DOUBLE PRECISION    :: eng,ens,pot
      REAL(KIND=GP), INTENT(IN)    :: dt
      INTEGER, INTENT(IN) :: t
      INTEGER             :: i,j,k

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
         WRITE(1,10) (t-1)*dt,eng,ens,pot
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
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a
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
! Computes the power spectrum
!
      tmp = 1./real(n,kind=GP)**6
      IF (ista.eq.1) THEN
         DO j = 1,n
            DO k = 1,n
               kmn = int(abs(ka(k))+1)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  Ek(kmn) = Ek(kmn)+tmp*abs(a(k,j,1))**2
               ENDIF
            END DO
         END DO
         DO i = 2,iend
            DO j = 1,n
               DO k = 1,n
                  kmn = int(abs(ka(k))+1)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     Ek(kmn) = Ek(kmn)+2*tmp*abs(a(k,j,i))**2
                  ENDIF
               END DO
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               DO k = 1,n
                  kmn = int(abs(ka(k))+1)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     Ek(kmn) = Ek(kmn)+2*tmp*abs(a(k,j,i))**2
                  ENDIF
               END DO
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
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b
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
      tmp = 1./real(n,kind=GP)**6
      IF (ista.eq.1) THEN
         DO j = 1,n
            DO k = 1,n
               kmn = int(sqrt(ka2(k,j,1))+.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  Ek(kmn) = Ek(kmn)+tmp*real(a(k,j,1)*conjg(b(k,j,1)))
               ENDIF
            END DO
         END DO
         DO i = 2,iend
            DO j = 1,n
               DO k = 1,n
                  kmn = int(sqrt(ka2(k,j,i))+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     Ek(kmn) = Ek(kmn)+2*tmp*real(a(k,j,i)*conjg(b(k,j,i)))
                  ENDIF
               END DO
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               DO k = 1,n
                  kmn = int(sqrt(ka2(k,j,i))+.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     Ek(kmn) = Ek(kmn)+2*tmp*real(a(k,j,i)*conjg(b(k,j,i)))
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
         OPEN(1,file='stransfer.' // nmb // '.txt')
         WRITE(1,30) Ektot
   30    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE sctrans
