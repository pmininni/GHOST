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
      
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend)  :: a,b
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,ista:iend) :: c
      COMPLEX(KIND=GP), DIMENSION(n,ista:iend) :: c1,c2
      REAL(KIND=GP), DIMENSION(n,jsta:jend)    :: r1,r2,rp
      REAL(KIND=GP) :: tmp
      INTEGER :: i,j

      tmp = 1.0_GP/real(n,kind=GP)**4

      DO j = jsta,jend
         DO i = 1,n
            c1(i,j) = a(i,j)
            c2(i,j) = b(i,j)
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
      
!**********************************************************************
