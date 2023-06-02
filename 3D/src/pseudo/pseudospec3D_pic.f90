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
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.ge.nth) private (k)
        DO j = 1,ny
          DO k = 1,nz
            phi(k,j,i) = rhoc(k,j,i)/(kk2(k,j,i)+kde2)
         END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE poisson_elecstat
