!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Extra subroutines for velocity formulation.
! You should use the FFTPLANS 
! and MPIVARS modules (see the file 'fftp2D_mod.f90') in each 
! program that call any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!
! 2012 Patricio Clark
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: patoclark@gmail.com
!=================================================================

!*****************************************************************
   SUBROUTINE elemwise2(a,b,c)
!-----------------------------------------------------------------
!
! Computes element-wise multiplication for two-dimensional
! matrices in real space.
!
! Parameters
!     a: input matrix
!     b: input matrix
!     c: element-wise product in Fourier space [output]

      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      IMPLICIT NONE
      
      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(n,ista:iend)  :: a,b
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,ista:iend) :: c
      COMPLEX(KIND=GP), DIMENSION(n,ista:iend) :: c1,c2
      REAL(KIND=GP), DIMENSION(n,jsta:jend)    :: r1,r2,rp
      REAL(KIND=GP) :: tmp
      INTEGER :: i,j

      tmp = 1.0_GP/real(n,kind=GP)**4

      DO i = ista,iend
         DO j = 1,n
            c1(j,i) = a(j,i)
            c2(j,i) = b(j,i)
         END DO
      END DO
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)

      DO j = jsta,jend
         DO i = 1,n
            rp(i,j) = r1(i,j)*r2(i,j)*tmp
         END DO
      END DO

      CALL fftp2d_real_to_complex(planrc,rp,c,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE elemwise2


!*****************************************************************
    SUBROUTINE gradre2(a,b,d,e)
!-----------------------------------------------------------------
!
! Two-dimensional inner product A.grad(A) in 
! real space. The components of the field A are 
! given by the matrixes a and b.
!
! Parameters
!     a: input matrix in the x-direction
!     b: input matrix in the y-direction
!     d: product (A.grad)A_x in Fourier space [output]
!     e: product (A.grad)A_y in Fourier space [output]
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend)  :: a,b
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,ista:iend) :: d,e
      COMPLEX(KIND=GP), DIMENSION(n,ista:iend) :: c1,c2,c3
      REAL(KIND=GP), DIMENSION(n,jsta:jend)    :: r1,r2,r3
      REAL(KIND=GP), DIMENSION(n,jsta:jend)    :: rx,ry
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j
      
      tmp = 1.0_GP/real(n,kind=GP)**4
!
! Computes (A_x.dx)A_dir
!
      c1 = a
      CALL derivk2(a,c2,1)
      CALL derivk2(b,c3,1)
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c3,r3,MPI_COMM_WORLD)

      DO j = jsta,jend
         DO i = 1,n
            rx(i,j) = r1(i,j)*r2(i,j)
            ry(i,j) = r1(i,j)*r3(i,j)
         END DO
      END DO
!
! Computes (A_y.dy)A_dir
!
      c1 = b
      CALL derivk2(a,c2,2)
      CALL derivk2(b,c3,2)
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c3,r3,MPI_COMM_WORLD)

      DO j = jsta,jend
         DO i = 1,n
            rx(i,j) = (rx(i,j)+r1(i,j)*r2(i,j))*tmp
            ry(i,j) = (ry(i,j)+r1(i,j)*r3(i,j))*tmp
         END DO
      END DO

      CALL fftp2d_real_to_complex(planrc,rx,d,MPI_COMM_WORLD)
      CALL fftp2d_real_to_complex(planrc,ry,e,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE gradre2
      
!*****************************************************************
      SUBROUTINE v_entrans(a,b,c,d,nmb)
!-----------------------------------------------------------------
!
! Computes the energy (or cross-helicity) transfer 
! in Fourier space in 2D. The output is written to 
! a file by the first node.
!
! Parameters
!     a  : field component in the x-direction
!     b  : field component in the y-direction
!     c  : nonlinear term in the x-direction
!     d  : nonlinear term in the y-direction
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
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: c,d
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
! Computes the energy flux
!
      tmp = 1.0_GP/real(n,kind=GP)**4
      IF (ista.eq.1) THEN
         DO j = 1,n
            kmn = int(sqrt(ka2(j,1))+.501)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
               Ek(kmn) = Ek(kmn)+                       &
                         (real(a(j,1)*conjg(c(j,1)))+   &
                          real(b(j,1)*conjg(d(j,1))))*tmp
            ENDIF
         END DO
         DO i = 2,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF (kmn.le.n/2+1) THEN
                  Ek(kmn) = Ek(kmn)+                         &
                            2*(real(a(j,i)*conjg(c(j,i)))+   &
                               real(b(j,i)*conjg(d(j,i))))*tmp
               ENDIF
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF (kmn.le.n/2+1) THEN
                  Ek(kmn) = Ek(kmn)+                         &
                            2*(real(a(j,i)*conjg(c(j,i)))+   &
                               real(b(j,i)*conjg(d(j,i))))*tmp
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
         OPEN(1,file='ktransfer.' // nmb // '.txt')
         WRITE(1,20) Ektot
   20    FORMAT( E23.15 )
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE v_entrans
      
!*****************************************************************
      SUBROUTINE vhdcheck(a,b,c,d,t)
!-----------------------------------------------------------------
!
! Computes
!
! Parameters
!     a  : velocity in the x-direction
!     b  : velocity in the y-direction
!     c  : force in the x-direction
!     d  : force in the y-direction
!     t  : time
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a,b,c,d
      DOUBLE PRECISION :: tmp,pot,u2,w2
      REAL(KIND=GP) :: t
      REAL(KIND=GP) :: tmq
      INTEGER       :: i,j

      tmq = 1.0_GP/real(n,kind=GP)**4
!
! Computes the mean energy and enstrophy
!
      CALL v_energy(a,b,u2,1)
      CALL v_energy(a,b,w2,0)
!
! Computes u.f
!
      tmp = 0.0D0
      IF (ista.eq.1) THEN
         DO j = 1,n
            tmp = tmp+(real(a(j,1)*conjg(c(j,1)))+   &
                       real(b(j,1)*conjg(d(j,1))))*tmq
         END DO
         DO i = 2,iend
            DO j = 1,n
               tmp = tmp+2*(real(a(j,i)*conjg(c(j,i)))+   &
                            real(b(j,i)*conjg(d(j,i))))*tmq
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               tmp = tmp+2*(real(a(j,i)*conjg(c(j,i)))+   &
                            real(b(j,i)*conjg(d(j,i))))*tmq
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
         WRITE(1,10) t,u2,w2,pot
   10    FORMAT( E13.6,E26.18,E26.18,E26.18 )
         CLOSE(1)
      ENDIF      

      RETURN
      END SUBROUTINE vhdcheck
      
!*****************************************************************
      SUBROUTINE v_energy(a,b,c,kin)
!-----------------------------------------------------------------
!
! Computes the mean energy or enstrophy in 2D.
! The output is valid only in the first node.
!
! Parameters
!     a  : velocity in the x-direction
!     b  : velocity in the y-direction
!     c  : at the output contains the energy or the enstrophy
!     kin: =1 computes the square of the velocity field
!          =0 computes the square of the vorticity_z
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a,b
      COMPLEX(KIND=GP), DIMENSION(n,ista:iend) :: c1,c2
      DOUBLE PRECISION, INTENT(OUT) :: c
      DOUBLE PRECISION              :: bloc
      REAL(KIND=GP)                 :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j

      bloc = 0.0D0
      tmp = 1.0_GP/real(n,kind=GP)**4
      
!
! Computes the square of the velocity field
!
      IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
            DO j = 1,n
               bloc = bloc+(abs(a(j,1))**2+abs(b(j,1))**2)*tmp
            END DO
            DO i = 2,iend
               DO j = 1,n
                  bloc = bloc+2*(abs(a(j,i))**2+abs(b(j,i))**2)*tmp
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  bloc = bloc+2*(abs(a(j,i))**2+abs(b(j,i))**2)*tmp
               END DO
            END DO
         ENDIF
!
! Computes the square of the vorticity_z
!
      ELSE IF (kin.eq.0) THEN
         CALL derivk2(a,c1,2)
         CALL derivk2(b,c2,1)
         DO i = ista,iend
            DO j = 1,n
               c1(j,i) = c2(j,i) - c1(j,i)
            END DO
         END DO
         
         IF (ista.eq.1) THEN
            DO j = 1,n
               bloc = bloc+abs(c1(j,1))**2*tmp
            END DO
            DO i = 2,iend
               DO j = 1,n
                  bloc = bloc+2*abs(c1(j,i))**2*tmp
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  bloc = bloc+2*abs(c1(j,i))**2*tmp
               END DO
            END DO
         ENDIF
      END IF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(bloc,c,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE v_energy
 
!*****************************************************************
      SUBROUTINE vspectrum(a,b,nmb)
!-----------------------------------------------------------------
!
! Computes the energy power spectrum in 2D. 
! The output is written to a file by the first node.
!
! Parameters
!     a  : input matrix in the x-directin
!     b  : input matrix in the y-direction
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

      tmp = 1.0_GP/real(n,kind=GP)**4
      
!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.0D0
      END DO  
!
! Computes the velocity square spectrum
!
      IF (ista.eq.1) THEN
         DO j = 1,n
            kmn = int(sqrt(ka2(j,1))+.5)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
               Ek(kmn) = Ek(kmn) + (abs(a(j,1))**2+abs(b(j,1))**2)*tmp
            ENDIF
         END DO
         DO i = 2,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  Ek(kmn) = Ek(kmn) + 2*(abs(a(j,i))**2+abs(b(j,i))**2)*tmp
               ENDIF
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  Ek(kmn) = Ek(kmn) + 2*(abs(a(j,i))**2+abs(b(j,i))**2)*tmp
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
      END SUBROUTINE vspectrum

