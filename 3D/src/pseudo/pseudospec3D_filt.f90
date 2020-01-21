!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Defines filters for different LES implementations. 
! You should use the FFTPLANS and MPIVARS modules (see the 
! file 'fftp_mod.f90') in each program that call any of the 
! subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2007 Jonathan Pietarila Graham and Pablo D. Mininni
!      National Center for Atmospheric Research.
!=================================================================

!*****************************************************************
      SUBROUTINE smooth3(a,b,c,d,e,f,alp)
!-----------------------------------------------------------------
!
! Alpha-smooths fields in three dimensions in Fourier
! space using a Helmholtz filter
!
! Parameters
!     a  : unsmoothed field in the x-direction
!     b  : unsmoothed field in the y-direction
!     c  : unsmoothed field in the z-direction
!     d  : smoothed field in the x-direction [output]
!     e  : smoothed field in the y-direction [output]
!     f  : smoothed field in the z-direction [output]
!     alp: value of alpha
!
      USE fprecision
      USE ali
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,nx,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,nx,ista:iend) :: d,e,f
      REAL(KIND=GP), INTENT(IN) :: alp
      REAL(KIND=GP)             :: tmp
      INTEGER                   :: i,j,k

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               IF ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) THEN
                  tmp = 1./(1.+alp**2*kk2(k,j,i))
                  d(k,j,i) = a(k,j,i)*tmp
                  e(k,j,i) = b(k,j,i)*tmp
                  f(k,j,i) = c(k,j,i)*tmp
               ELSE
                  d(k,j,i) = 0.0_GP
                  e(k,j,i) = 0.0_GP
                  f(k,j,i) = 0.0_GP
               ENDIF
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE smooth3

!*****************************************************************
      SUBROUTINE smooth(a,b,alp)
!-----------------------------------------------------------------
!
! Alpha-smooths one scalar field in Fourier space using a
! Helmholtz filter
!
! Parameters
!     a  : unsmoothed field
!     b  : smoothed field [output]
!     alp: value of alpha
!
      USE fprecision
      USE ali
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: b
      REAL(KIND=GP), INTENT(IN) :: alp
      REAL(KIND=GP)             :: tmp
      INTEGER                   :: i,j,k

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               IF ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) THEN
                  tmp = 1./(1.+alp**2*kk2(k,j,i))
                  b(k,j,i) = a(k,j,i)*tmp
               ELSE
                  b(k,j,i) = 0.0_GP
               ENDIF
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE smooth
      
!*****************************************************************
      SUBROUTINE gaussian(a,b,delta)
!-----------------------------------------------------------------
!
! Gaussian filter in Fourier space. The characteristic cut-off
! wavenumber is k_c = pi/Delta (see S. Pope, "Turbulent Flows,"
! Chap. 13.2).
!
! Parameters
!     a    : unsmoothed field
!     b    : smoothed field [output]
!     delta: Cut-off length
!
      USE fprecision
      USE ali
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: b
      REAL(KIND=GP), INTENT(IN) :: delta
      REAL(KIND=GP)             :: tmp
      INTEGER                   :: i,j,k

      tmp = delta**2/24
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               IF (kn2(k,j,i).le.kmax) THEN
                  b(k,j,i) = a(k,j,i)*exp(-kk2(k,j,i)*tmp)
               ELSE
                  b(k,j,i) = 0.0_GP
               ENDIF
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE gaussian

!*****************************************************************
      SUBROUTINE box(a,b,delta)
!-----------------------------------------------------------------
!
! Box filter in Fourier space. The characteristic cut-off
! wavenumber is k_c = pi/Delta (see S. Pope, "Turbulent Flows,"
! Chap. 13.2).
!
! Parameters
!     a    : unsmoothed field
!     b    : smoothed field [output]
!     delta: Cut-off length
!
      USE fprecision
      USE ali
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: b
      REAL(KIND=GP), INTENT(IN) :: delta
      REAL(KIND=GP)             :: tmp
      INTEGER                   :: i,j,k

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               IF (kn2(k,j,i).le.kmax) THEN
                  tmp = .5_GP*delta*sqrt(kk2(k,j,i))
                  b(k,j,i) = a(k,j,i)*sin(tmp)/tmp
               ELSE
                  b(k,j,i) = 0.0_GP
               ENDIF
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE box

!*****************************************************************
      SUBROUTINE sharp(a,b,delta)
!-----------------------------------------------------------------
!
! Sharp filter in Fourier space. The characteristic cut-off
! wavenumber is k_c = pi/Delta (see S. Pope, "Turbulent Flows,"
! Chap. 13.2).
!
! Parameters
!     a    : unsmoothed field
!     b    : smoothed field [output]
!     delta: Cut-off length
!
      USE fprecision
      USE var
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: b
      REAL(KIND=GP), INTENT(IN) :: delta
      REAL(KIND=GP)             :: tmp
      INTEGER                   :: i,j,k

      tmp = (pi/delta)**2
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               IF (kn2(k,j,i).le.tmp) THEN
                  b(k,j,i) = a(k,j,i)
               ELSE
                  b(k,j,i) = 0.0_GP
               ENDIF
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE sharp
