!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines to compute spatial derivatives and nonlinear 
! terms in surface quasi-geostrophic (SQG) equations using
! a pseudo-spectral method. You should use the FFTPLANS and
! MPIVARS modules (see the file 'fftp2D_mod.f90') in each 
! program that call any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!
! 2011 Tomas Teitelbaum.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: teitelbaum@df.uba.ar 
!=================================================================

!*****************************************************************
      SUBROUTINE scalarder(a,b)
!-----------------------------------------------------------------
!
! Computes the scalar stream function in SQG
!
! Parameters
!     a: input matrix
!     b: the output is the scalar function q
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
            b(j,i) = -sqrt(ka2(j,i))*a(j,i)
         END DO
      END DO

      RETURN
      END SUBROUTINE scalarder


!*****************************************************************
      SUBROUTINE sqgenergy(a,b,kin)
!-----------------------------------------------------------------
!
! Computes the mean kinetic energy, the variance of 
! the scalar field, and the square vorticity. The 
! output is valid only in the first node.
!
! Parameters
!     a  : input matrix with the scalar field
!     b  : the output contains the energy
!     kin: =2 computes the square of the scalar field
!          =1 computes the energy
!          =0 computes the vorticity
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
               bloc = bloc+sqrt(ka2(j,1))*abs(a(j,1))**2*tmp
            END DO
            DO i = 2,iend
               DO j = 1,n
                  bloc = bloc+2*sqrt(ka2(j,i))*abs(a(j,i))**2*tmp
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  bloc = bloc+2*sqrt(ka2(j,i))*abs(a(j,i))**2*tmp
               END DO
            END DO
         ENDIF
!
! Computes the vorticity
!
      ELSE
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
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(bloc,b,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE sqgenergy

!*****************************************************************
      SUBROUTINE sqgcheck(a,b,t)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of energy in SQG
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
      CALL sqgenergy(a,eng,1)
      CALL sqgenergy(a,ens,0)
!
! Computes the energy injection rate
!
      tmp = 0.0D0
      IF (ista.eq.1) THEN
         DO j = 1,n
            tmp = tmp+sqrt(ka2(j,1))*real(b(j,1)*conjg(a(j,1)))*tmq
         END DO
         DO i = 2,iend
            DO j = 1,n
               tmp = tmp+2*sqrt(ka2(j,i))*real(b(j,i)*conjg(a(j,i)))*tmq
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               tmp = tmp+2*sqrt(ka2(j,i))*real(b(j,i)*conjg(a(j,i)))*tmq
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
      END SUBROUTINE sqgcheck

!*****************************************************************
      SUBROUTINE sqgspectrum(a,nmb)
!-----------------------------------------------------------------
!
! Computes the energy power spectrum in 2D. 
! The output is written to a file by the first node.
!
! Parameters
!     a  : streamfunction or vector potential
!     nmb: the extension used when writting the file
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
               Ek(kmn) = Ek(kmn)+sqrt(ka2(j,1))*abs(a(j,1))**2*tmp
            ENDIF
         END DO
         DO i = 2,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  Ek(kmn) = Ek(kmn)+2*sqrt(ka2(j,i))*abs(a(j,i))**2*tmp
               ENDIF
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  Ek(kmn) = Ek(kmn)+2*sqrt(ka2(j,i))*abs(a(j,i))**2*tmp
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
         OPEN(1,file='kspectrum.' // nmb // '.txt')
         WRITE(1,20) Ektot
   20    FORMAT( E23.15 )
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE sqgspectrum
