!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines for computing spatial derivatives and nonlinear 
! terms in incompressible MHD and Hall-MHD equations in 2D 
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
      SUBROUTINE poissonb0(a,b,c,b0)
!-----------------------------------------------------------------
!
! Poisson bracket of the scalar fields A and B 
! in real space.
!
! Parameters
!     a : input matrix
!     b : input matrix
!     c : Poisson bracket {A,B} [output]
!     b0: amplitude of the uniform field in y
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
      REAL(KIND=GP), INTENT(IN) :: b0
      REAL(KIND=GP) :: tmp
      INTEGER       :: i,j

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
! Computes (dA/dy+b0).dB/dx
!
      CALL derivk2(a,c1,2)
      CALL derivk2(b,c2,1)
      c2(1,1) = -b0*real(n,KIND=GP)**2
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
      END SUBROUTINE poissonb0

!*****************************************************************
      SUBROUTINE mhdcheck(a,b,c,d,t)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of energy in MHD 2D
!
! Parameters
!     a  : streamfunction
!     b  : vector potential
!     c  : external kinetic force
!     d  : external magnetic force
!     t  : time
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a,b,c,d
      DOUBLE PRECISION :: engk,engm,eng,asq,udb
      DOUBLE PRECISION :: potk,potm,ens,cur,tmp
      REAL(KIND=GP) :: t
      REAL(KIND=GP) :: tmq
      INTEGER       :: i,j

      tmq = 1.0_GP/real(n,kind=GP)**4
!
! Computes the mean energy, enstrophy, square 
! current, and square vector potential.
!
      CALL energy(a,engk,1)
      CALL energy(b,engm,1)
      CALL energy(a,ens,0)
      CALL energy(b,cur,0)
      CALL energy(b,asq,2)
      eng = engk+engm
!
! Computes the kinetic energy injection rate
!
      tmp = 0.0D0
      IF (ista.eq.1) THEN
         DO j = 1,n
            tmp = tmp+ka2(j,1)*real(c(j,1)*conjg(a(j,1)))*tmq
         END DO
         DO i = 2,iend
            DO j = 1,n
               tmp = tmp+2*ka2(j,i)*real(c(j,i)*conjg(a(j,i)))*tmq
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               tmp = tmp+2*ka2(j,i)*real(c(j,i)*conjg(a(j,i)))*tmq
            END DO
         END DO
      ENDIF
      CALL MPI_REDUCE(tmp,potk,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Computes the magnetic energy injection rate
!
      tmp = 0.0D0
      IF (ista.eq.1) THEN
         DO j = 1,n
            tmp = tmp+ka2(j,1)*real(d(j,1)*conjg(b(j,1)))*tmq
         END DO
         DO i = 2,iend
            DO j = 1,n
               tmp = tmp+2*ka2(j,i)*real(d(j,i)*conjg(b(j,i)))*tmq
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               tmp = tmp+2*ka2(j,i)*real(d(j,i)*conjg(b(j,i)))*tmq
            END DO
         END DO
      ENDIF
      CALL MPI_REDUCE(tmp,potm,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Computes the cross correlation between
! velocity and magnetic fields
!
      tmp = 0.
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
      CALL MPI_REDUCE(tmp,udb,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='balance.txt',position='append')
         WRITE(1,10) t,eng,ens,cur
   10    FORMAT( E13.6,E22.14,E22.14,E22.14 )
         CLOSE(1)
         OPEN(1,file='cross.txt',position='append')
         WRITE(1,10) t,udb,asq
   20    FORMAT( E13.6,E22.14,E22.14 )
         CLOSE(1)
         OPEN(1,file='energy.txt',position='append')
         WRITE(1,20) t,engk,engm
         CLOSE(1)
         OPEN(1,file='inject.txt',position='append')
         WRITE(1,20) t,potk,potm
         CLOSE(1)
      ENDIF      

      RETURN
      END SUBROUTINE mhdcheck

!*****************************************************************
      SUBROUTINE pmspectrum(ps,a,nmb,kin)
!-----------------------------------------------------------------
!
! Computes the energy power spectrum in 2D for E+/-. 
! The output is written to a file by the first node.
!
! Parameters
!     ps : streamfunction 
!     a  : vector potential
!     nmb: the extension used when writting the file
!     kin: =0 computes the E+ spectrum
!          =1 computes the E- spectrum
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: ps,a
      DOUBLE PRECISION, DIMENSION(n/2+1)  :: Ek,Ektot
      REAL(KIND=GP)       :: q, sgn, tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: kmn
      INTEGER             :: i,j
      CHARACTER(len=*), INTENT(IN) :: nmb

      sgn = 2.0_GP
      IF (kin.eq.1) THEN
        sgn = -2.0_GP
      ENDIF
!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.0_GP
      END DO
!
! Computes the energy spectrum
!
      tmp = 1.0_GP/real(n,kind=GP)**4
      IF (ista.eq.1) THEN
         DO j = 1,n
            kmn = int(sqrt(ka2(j,1))+.5)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
               q = abs(ps(j,1))**2+abs(a(j,1))**2+     &
                   sgn*real(ps(j,1)*conjg(a(j,1)))
               Ek(kmn) = Ek(kmn)+ka2(j,1)*q*tmp
            ENDIF
         END DO
         DO i = 2,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  q = abs(ps(j,i))**2+abs(a(j,i))**2+  & 
                      sgn*real(ps(j,i)*conjg(a(j,i)))
                  Ek(kmn) = Ek(kmn)+2*ka2(j,i)*q*tmp
               ENDIF
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  q = abs(ps(j,i))**2+abs(a(j,i))**2+  &
                      sgn*real(ps(j,i)*conjg(a(j,i)))
                  Ek(kmn) = Ek(kmn)+2*ka2(j,i)*q*tmp
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
            OPEN(1,file='epspectrum.' // nmb // '.txt')
         ELSE 
            OPEN(1,file='emspectrum.' // nmb // '.txt')
         ENDIF
         WRITE(1,20) Ektot
   20    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE pmspectrum

!*****************************************************************
      SUBROUTINE vectrans(a,b,nmb)
!-----------------------------------------------------------------
!
! Computes the square vector potential transfer in 
! Fourier space in 2D MHD. The output is written 
! to a file by the first node.
!
! Parameters
!     a  : vector potential
!     b  : Poisson bracket of the streamfunction and a
!     nmb: the extension used when writting the file
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE filefmt
      USE grid
      USE kes
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a,b
      DOUBLE PRECISION, DIMENSION(n/2+1) :: Ek,Ektot
      REAL(KIND=GP)       :: tmp
      INTEGER             :: kmn
      INTEGER             :: i,j
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.0_GP
      END DO
!
! Computes the square vector potential flux
!
      tmp = 1.0_GP/real(n,kind=GP)**4
      IF (ista.eq.1) THEN
         DO j = 1,n
            kmn = int(sqrt(ka2(j,1))+.501)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
              Ek(kmn) = Ek(kmn)+real(b(j,1)*conjg(a(j,1)))*tmp
            ENDIF
         END DO
         DO i = 2,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF (kmn.le.n/2+1) THEN
                  Ek(kmn) = Ek(kmn)+2*real(b(j,i)*conjg(a(j,i)))*tmp
               ENDIF
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF (kmn.le.n/2+1) THEN
                  Ek(kmn) = Ek(kmn)+2*real(b(j,i)*conjg(a(j,i)))*tmp
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
         OPEN(1,file='vectransf.' // nmb // '.txt')
         WRITE(1,20) Ektot
   20    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE vectrans
