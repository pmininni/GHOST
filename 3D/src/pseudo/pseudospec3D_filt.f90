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
!     e  : smoothed field in the x-direction [output]
!     f  : smoothed field in the x-direction [output]
!     alp: value of alpha
!
      USE ali
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX, INTENT(IN), DIMENSION(n,n,ista:iend)  :: a,b,c
      COMPLEX, INTENT(OUT), DIMENSION(n,n,ista:iend) :: d,e,f
      REAL, INTENT(IN) :: alp
      REAL             :: tmp
      INTEGER          :: i,j,k

      DO i = ista,iend
         DO j = 1,n
            DO k = 1,n
               IF ((ka2(k,j,i).le.kmax).and.(ka2(k,j,i).ge.tiny)) THEN
                  tmp = 1./(1.+alp**2*ka2(k,j,i))
                  d(k,j,i) = a(k,j,i)*tmp
                  e(k,j,i) = b(k,j,i)*tmp
                  f(k,j,i) = c(k,j,i)*tmp
               ELSE
                  d(k,j,i) = 0.
                  e(k,j,i) = 0.
                  f(k,j,i) = 0.
               ENDIF
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE smooth3
