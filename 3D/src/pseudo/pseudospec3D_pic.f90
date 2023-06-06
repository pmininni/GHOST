!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Extra subroutines to compute spatial derivatives and 
! nonlinear terms in Hybrid-PIC solvers using a
! pseudo-spectral method. You should use the FFTPLANS and 
! MPIVARS modules (see the file 'fftp_mod.f90') in each 
! program that calls any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2023 Facundo Pugliese
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: fpugliese@df.uba.ar 
!=================================================================

!*****************************************************************
      SUBROUTINE poisson_elecstat(rhoc,kde2,phi)
!-----------------------------------------------------------------
! 
! Solves poisson problem lap(phi) = -rhoc + kde2*phi wherw kde2 is
! the square inverse Debye length.
!
! Parameters
!     rhoc : input complex density
!     kde2 : square inverse Debye length
!     phi  : output complex potential
!
      USE fprecision
      USE kes
      USE grid
      USE mpivars
!$    USE threads

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: rhoc
      REAL(KIND=GP)   , INTENT (IN)                             :: kde2
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: phi
      IF (ista.eq.0) THEN
         phi(1,1,1) = 0.0_GP
         DO k = 2,nz
            phi(k,1,1) = rhoc(k,1,1)/(kk2(k,1,1)+kde2)
         END DO
!$omp parallel do if (iend-ista.ge.nth) private (k)
         DO j = 2,ny
            DO k = 1,nz
               phi(k,j,1) = rhoc(k,j,1)/(kk2(k,j,1)+kde2)
            END DO
         END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = 2,iend
!$omp parallel do if (iend-ista.ge.nth) private (k)
           DO j = 1,ny
             DO k = 1,nz
               phi(k,j,i) = rhoc(k,j,i)/(kk2(k,j,i)+kde2)
            END DO
           END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.ge.nth) private (k)
           DO j = 1,ny
             DO k = 1,nz
               phi(k,j,i) = rhoc(k,j,i)/(kk2(k,j,i)+kde2)
            END DO
           END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE poisson_elecstat

!*****************************************************************
      SUBROUTINE ehpiccheck(Tem, Ex, Ey, Ez, t, dt)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of the total
! energy of electrostatic hybrid-pic solver.
!
! Output files contain:
! 'energy.txt':   time, <v^2>, <E^2>
!
! Parameters
!     Tem: kinetic energy (of particles) in real space
!     Ex : electric field in the x-direction
!     Ey : electric field in the y-direction
!     Ez : electric field in the z-direction
!     t  : number of time steps made
!     dt : time step
!
      USE fprecision
      USE commtypes
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: Ex,Ey,Ez
      REAL(KIND=GP)   , INTENT(IN), DIMENSION(ksta:kend,ny,nz) :: Tem
      REAL(KIND=GP), INTENT(IN)    :: dt
      INTEGER      , INTENT(IN)    :: t
      DOUBLE PRECISION    :: ek,ekloc,ep
      INTEGER             :: i,j,k
      REAL(KIND=GP)       :: rmp

      CALL energy(Ex,Ey,Ez,ep,1)
      
      rmp = 1/(real(nx,KIND=GP)*real(ny,KIND=GP)*real(nz,KIND=GP))
      ekloc = 0.0D0
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               ekloc = ekloc + Tem(k,j,i)*rmp
            ENDDO
         ENDDO
      ENDDO

      CALL MPI_REDUCE(ekloc,ek,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='energy.txt',position='append')
         WRITE(1,10) (t-1)*dt,ek,ep
   10    FORMAT( E13.6,E22.14,E22.14,E22.14 )
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE ehpiccheck

