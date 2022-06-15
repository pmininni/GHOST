!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines to compute nudging term and associated quantities.  The
! nudging solver overrides the usual forcing by using initialfv to read
! the reference files and then uses the option rand=2 to interpolate in
! time.  You should use the FFTPLANS and MPIVARS modules (see the file
! 'fftp_mod.f90') in each program that calls any of the subroutines in
! this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!      
! 2022 Patricio Clark Di Leoni
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: clark@df.uba.ar 
!=================================================================

!*****************************************************************
      SUBROUTINE ndg_term(a,b,c,amp,kdn,kup)
!-----------------------------------------------------------------
!
! Computes the nudging term -amp*(a-b).
! A band pass filter is applied between kdn and kup.
!
! Parameters
!     a  : first field
!     b  : second field
!     c  : at the output contains the nudging term
!     amp: amplitude of the nudging term
!     kdn: minimum wave number
!     kup: maximum wave number
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend)  :: a,b
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: c
      REAL(KIND=GP), INTENT(IN) :: amp, kdn, kup
      INTEGER             :: i,j,k

      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
      DO j = 1,ny
      DO k = 1,nz
         IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
            c(k,j,i) = -amp*(a(k,j,i) - b(k,j,i))
         ELSE
            c(k,j,i) = 0.0_GP
         ENDIF
      END DO
      END DO
      END DO

      RETURN
      END SUBROUTINE ndg_term
