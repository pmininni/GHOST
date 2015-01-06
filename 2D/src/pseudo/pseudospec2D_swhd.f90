!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Extra subroutines stuff in shallow waters.
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
      SUBROUTINE swhdcheck(a,b,c,d,e,f,g,t)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of energy in SWHD
!
! Parameters
!     a  : velocity in the x-direction
!     b  : velocity in the y-direction
!     c  : height
!     d  : force in the x-direction
!     e  : force in the y-direction
!     f  : bottom topography
!     g  : gravity
!     t  : time
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), DIMENSION(n,ista:iend) :: c1,c2,c3,c4
      DOUBLE PRECISION :: eng,eps,zet,cin,pot,tmp,tmr
      DOUBLE PRECISION :: var,grd
      REAL(KIND=GP) :: g,t
      REAL(KIND=GP) :: tmq
      INTEGER       :: i,j

      tmq = 1.0_GP/real(n,kind=GP)**4

!
! Computes kinetic energy and the gravitational potential energy
!
      CALL swenergy(a,b,c,f,cin,g,1)
      CALL swenergy(a,b,c,f,pot,g,2)
!
! Computes the total energy
!
      eng = cin + pot
!
! Computes the energy injection rate
!
      tmp = 0.0D0
      CALL elemwise2(a,c,c1)
      CALL elemwise2(b,c,c2)
      IF (ista.eq.1) THEN
         DO j = 1,n
            tmp = tmp+(real(c1(j,1)*conjg(d(j,1)))+   &
                       real(c2(j,1)*conjg(e(j,1))))*tmq
         END DO
         DO i = 2,iend
            DO j = 1,n
               tmp = tmp+2*(real(c1(j,i)*conjg(d(j,i)))+   &
                            real(c2(j,i)*conjg(e(j,i))))*tmq
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               tmp = tmp+2*(real(c1(j,i)*conjg(d(j,i)))+   &
                            real(c2(j,i)*conjg(e(j,i))))*tmq
            END DO
         END DO
      ENDIF
      CALL MPI_REDUCE(tmp,eps,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Computes the energy dissipation rate
!
      tmr = 0.0D0
      CALL laplak2(a,c3)
      CALL laplak2(b,c4)
      IF (ista.eq.1) THEN
         DO j = 1,n
            tmr = tmr+(real(c1(j,1)*conjg(c3(j,1)))+   &
                       real(c2(j,1)*conjg(c4(j,1))))*tmq
         END DO
         DO i = 2,iend
            DO j = 1,n
               tmr = tmr+2*(real(c1(j,i)*conjg(c3(j,i)))+   &
                            real(c2(j,i)*conjg(c4(j,i))))*tmq
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               tmr = tmr+2*(real(c1(j,i)*conjg(c3(j,i)))+   &
                            real(c2(j,i)*conjg(c4(j,i))))*tmq
            END DO
         END DO
      ENDIF
      CALL MPI_REDUCE(tmr,zet,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Computes the variance and the gradient of the scalar
!
      CALL variance(c,var,1)
      CALL variance(c,grd,0)
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='swbalance.txt',position='append')
         WRITE(1,10) t,eng,zet,eps
   10    FORMAT( E13.6,E26.18,E26.18,E26.18 )
         CLOSE(1)
         OPEN(1,file='energy.txt',position='append')
         WRITE(1,20) t,cin,pot
   20    FORMAT( E13.6,E26.18,E26.18 )
         CLOSE(1)
         OPEN(1,file='scalar.txt',position='append')
         WRITE(1,30) t,var,grd
   30    FORMAT( E13.6,E22.14,E22.14 )
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE swhdcheck

!*****************************************************************
      SUBROUTINE swenergy(a,b,c,d,e,g,kin)
!-----------------------------------------------------------------
!
! Computes the mean energy or enstrophy in 2D.
! The output is valid only in the first node.
!
! Parameters
!     a  : velocity in the x-direction
!     b  : velocity in the y-direction
!     c  : height
!     d  : bottom topography
!     e  : output
!     g  : gravity
!     kin: =1 computes the kinetic energy
!          =2 computes the gravitational potential energy
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      USE kes
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a,b,c,d
      REAL(KIND=GP), DIMENSION(n,jsta:jend) :: r1,r2,rh
      COMPLEX(KIND=GP), DIMENSION(n,ista:iend) :: cb,c1,c2
      DOUBLE PRECISION, INTENT(OUT) :: e
      DOUBLE PRECISION              :: bloc
      REAL(KIND=GP)                 :: g
      REAL(KIND=GP)                 :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j

      bloc = 0.0D0
      tmp = 1.0_GP/real(n,kind=GP)**2
      
!
! Computes the kinetic energy
!
      IF (kin.eq.1) THEN
         DO i = ista,iend
            DO j = 1,n
               cb(j,i) = c(j,i) - d(j,i)
            END DO
         END DO
         CALL swenergyaux(a,cb,r1)
         CALL swenergyaux(b,cb,r2)
         
         DO j = jsta,jend
            DO i = 1,n
               bloc = bloc + (r1(i,j) + r2(i,j))*tmp
            END DO
         END DO
!
! Computes the gravitational potential energy
!
      ELSE IF (kin.eq.2) THEN
         IF (ista.eq.1) THEN
            DO j = 1,n
               bloc = bloc+g*abs(c(j,1))**2*tmp**2
            END DO
            DO i = 2,iend
               DO j = 1,n
                  bloc = bloc+2*g*abs(c(j,i))**2*tmp**2
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  bloc = bloc+2*g*abs(c(j,i))**2*tmp**2
               END DO
            END DO
         ENDIF
      END IF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(bloc,e,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE swenergy
      
!*****************************************************************
      SUBROUTINE swenergyaux(a,b,c)
!-----------------------------------------------------------------
!
! Computes auxiliary variable used to calculate the energy
! and the enstrophy
!
! Parameters
!     a  : input matrix
!     b  : input matrix
!     c  : output
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend)  :: a,b
      REAL(KIND=GP), INTENT(OUT), DIMENSION(n,jsta:jend) :: c
      COMPLEX(KIND=GP), DIMENSION(n,ista:iend) :: c1,c2
      REAL(KIND=GP), DIMENSION(n,jsta:jend)    :: r1,r2
      REAL(KIND=GP)       :: tmp
      INTEGER             :: i,j

      tmp = 1.0_GP/real(n,kind=GP)**2
      
      DO i = ista,iend
         DO j = 1,n
            c1(j,i) = a(j,i)*tmp
            c2(j,i) = b(j,i)*tmp
         END DO
      END DO
      
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
            
      DO j = jsta,jend
         DO i = 1,n
            c(i,j) = (r1(i,j)**2)*r2(i,j)
         END DO
      END DO

      RETURN
      END SUBROUTINE swenergyaux      

!*****************************************************************
      SUBROUTINE swspectrum(a,b,c,d,g,nmb,kin)
!-----------------------------------------------------------------
!
! Computes the energy power spectrum in 2D. 
! The output is written to a file by the first node.
!
! Parameters
!     a  : input matrix in the x-directin
!     b  : input matrix in the y-direction
!     c  : height
!     d  : bottom topography
!     g  : gravity
!     nmb: the extension used when writting the file
!     kin: =1 computes the energy spectrum
!          =0 computes the enstrophy spectrum
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE filefmt
      USE grid
      USE kes
      USE ali
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a,b,c,d
      COMPLEX(KIND=GP), DIMENSION(n,ista:iend) :: cb,c1,c2
      DOUBLE PRECISION, DIMENSION(n/2+1) :: Ek,Ektot
      REAL(KIND=GP)       :: g
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: kmn
      INTEGER             :: i,j
      CHARACTER(len=*), INTENT(IN) :: nmb

      tmp = 1.0_GP/real(n,kind=GP)**4
      
      cb = c - d

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.0D0
      END DO  
! 
! Computes the kinetic energy spectrum
! 
      IF (kin.eq.1) THEN
         CALL swspectrumaux(a,cb,c1)
         CALL swspectrumaux(b,cb,c2)
         
         IF (ista.eq.1) THEN
            DO j = 1,n
               kmn = int(sqrt(ka2(j,1))+.5)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  Ek(kmn) = Ek(kmn)+(abs(c1(j,1))**2+abs(c2(j,1))**2)*tmp
               ENDIF
            END DO
            DO i = 2,iend
               DO j = 1,n
                  kmn = int(sqrt(ka2(j,i))+.5)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     Ek(kmn) = Ek(kmn)+2*(abs(c1(j,i))**2+abs(c2(j,i))**2)*tmp
                  ENDIF
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  kmn = int(sqrt(ka2(j,i))+.5)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     Ek(kmn) = Ek(kmn)+2*(abs(c1(j,i))**2+abs(c2(j,i))**2)*tmp
                  ENDIF
               END DO
            END DO
         ENDIF
! 
! Computes the gravitational potential energy spectrum
! 
      ELSE IF (kin.eq.2) THEN
         IF (ista.eq.1) THEN
            DO j = 1,n
               kmn = int(sqrt(ka2(j,1))+.5)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  Ek(kmn) = Ek(kmn) + g*abs(c(j,1))**2*tmp
               ENDIF
            END DO
            DO i = 2,iend
               DO j = 1,n
                  kmn = int(sqrt(ka2(j,i))+.5)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     Ek(kmn) = Ek(kmn) + 2*g*abs(c(j,i))**2*tmp
                  ENDIF
               END DO
            END DO
         ELSE
            DO i = ista,iend
               DO j = 1,n
                  kmn = int(sqrt(ka2(j,i))+.5)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     Ek(kmn) = Ek(kmn) + 2*g*abs(c(j,i))**2*tmp
                  ENDIF
               END DO
            END DO
         ENDIF
      END IF

!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF (kin.eq.1) THEN
            OPEN(1,file='kinspectrum.' // nmb // '.txt')
         ELSE IF (kin.eq.2) THEN
            OPEN(1,file='graspectrum.' // nmb // '.txt')
         ENDIF
         WRITE(1,20) Ektot
   20    FORMAT( E23.15 )
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE swspectrum
      

!*****************************************************************
      SUBROUTINE swspectrumaux(a,b,c)
!-----------------------------------------------------------------
!
! Computes auxiliary variable used to calculate the kinetic energy
! and the gravitational potential energy spectrum
!
! Parameters
!     a  : input matrix
!     b  : input matrix
!     c  : output
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      USE kes
      USE ali
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend)  :: a,b
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,ista:iend) :: c
      COMPLEX(KIND=GP), DIMENSION(n,ista:iend) :: c1,c2
      REAL(KIND=GP), DIMENSION(n,jsta:jend)    :: r1,r2,r3
      REAL(KIND=GP)       :: tmp
      INTEGER             :: i,j,i2,j2

      tmp = 1.0_GP/real(n,kind=GP)**2
      
      DO i = ista,iend
         DO j = 1,n
            c1(j,i) = a(j,i)*tmp
            c2(j,i) = b(j,i)*tmp
         END DO
      END DO
      
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      
      DO j = jsta,jend
         DO i = 1,n
            r3(i,j) = SQRT((r1(i,j)**2)*r2(i,j))
         END DO
      END DO

      CALL fftp2d_real_to_complex(planrc,r3,c,MPI_COMM_WORLD)

      ! De-aliasing
      DO i = ista,iend
         DO j = 1,n
            IF (ka2(j,i).gt.kmax) c(j,i) = 0
         END DO
      END DO
      
      RETURN
      END SUBROUTINE swspectrumaux
     
!*****************************************************************
      SUBROUTINE point(a,t)
!-----------------------------------------------------------------
!
! Writes the value of the field in the (1,1) position
!
! Parameters:
!     a  : input matrix
!     t  : time
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend)  :: a
      COMPLEX(KIND=GP), DIMENSION(n,ista:iend) :: c
      REAL(KIND=GP),    DIMENSION(n,jsta:jend) :: r
      REAL(KIND=GP)       :: t,tmp
      INTEGER             :: i,j

      tmp = 1.0_GP/real(n,kind=GP)**2
      
      DO i = ista,iend
         DO j = 1,n
            c(j,i) = a(j,i)*tmp
         END DO
      END DO
      
      CALL fftp2d_complex_to_real(plancr,c,r,MPI_COMM_WORLD)
            
      IF (myrank.eq.0) THEN
         OPEN(1,file='point.txt',position='append')
         WRITE(1,10) t,r(1,1)
   10    FORMAT( E13.6,E26.18 )
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE point

!*****************************************************************
      SUBROUTINE buhler(a,b,c,d,e)
!-----------------------------------------------------------------
!
!     Calculates grad(h)*grad(A)/h
!
! Parameters
!     a: input matrix
!     b: input matrix in the x direction
!     c: input matrix in the y direction
!     d: (grad(a)*grad(b)/a)_x in Fourier space
!     e: (grad(a)*grad(b)/a)_y in Fourier space
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend)  :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,ista:iend) :: d,e
      COMPLEX(KIND=GP), DIMENSION(n,ista:iend) :: c1,c2,c3,c4
      REAL(KIND=GP), DIMENSION(n,jsta:jend)    :: r1,r2,r3,r4
      REAL(KIND=GP), DIMENSION(n,jsta:jend)    :: rx,ry
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j
      
      tmp = 1.0_GP/real(n,kind=GP)**2
!
!     (grad(h)*grad(A)/h)_x
!
      c1 = a
      CALL derivk2(a,c2,1)
      CALL derivk2(b,c3,1)
      CALL derivk2(c,c4,1)
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c3,r3,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c4,r4,MPI_COMM_WORLD)

      DO j = jsta,jend
         DO i = 1,n
            rx(i,j) = r2(i,j)*r3(i,j)
            ry(i,j) = r2(i,j)*r4(i,j)
         END DO
      END DO
!
!     (grad(h)*grad(A)/h)_y
!
      CALL derivk2(a,c2,2)
      CALL derivk2(b,c3,2)
      CALL derivk2(c,c4,2)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c3,r3,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c4,r4,MPI_COMM_WORLD)

      DO j = jsta,jend
         DO i = 1,n
            rx(i,j) = (rx(i,j)+r2(i,j)*r3(i,j))*tmp/r1(i,j)
            ry(i,j) = (ry(i,j)+r2(i,j)*r4(i,j))*tmp/r1(i,j)
         END DO
      END DO

      CALL fftp2d_real_to_complex(planrc,rx,d,MPI_COMM_WORLD)
      CALL fftp2d_real_to_complex(planrc,ry,e,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE buhler
