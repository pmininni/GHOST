!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Extra subroutines to compute nonlinear terms in the 
! incompressible alpha-model HD, MHD and Hall-MHD equations 
! in 3D using a pseudo-spectral method. 
! You should use the FFTPLANS and MPIVARS modules (see the 
! file 'fftp_mod.f90') in each program that call any of the 
! subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2004 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.edu
!=================================================================

!*****************************************************************
      SUBROUTINE aprodre3(a,b,c,d,e,f,alp)
!-----------------------------------------------------------------
!
! Three-dimensional cross product curl(A)xAs in 
! real space. The components of the field A are 
! given by the matrixes a,b and c
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     d  : product [curl(A)xAs]_x in Fourier space [output]
!     e  : product [curl(A)xAs]_y in Fourier space [output]
!     f  : product [curl(A)xAs]_z in Fourier space [output]
!     alp: value of alpha
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend)  :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,n,ista:iend) :: d,e,f
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)    :: r1,r2
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)    :: r3,r4
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)    :: r5,r6
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)    :: r7
      REAL(KIND=GP), INTENT(IN) :: alp
      REAL(KIND=GP)             :: tmp
      INTEGER          :: i,j,k

!
! Computes curl(A)
!
      CALL rotor3(b,c,d,1)
      CALL rotor3(a,c,e,2)
      CALL rotor3(a,b,f,3)
      CALL fftp3d_complex_to_real(plancr,d,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,e,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,f,r3,MPI_COMM_WORLD)
!
! Computes As
!
      CALL smooth3(a,b,c,d,e,f,alp)
      CALL fftp3d_complex_to_real(plancr,d,r4,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,e,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,f,r6,MPI_COMM_WORLD)
!
! Computes curl(A)xAs
!
      tmp = 1./real(n,kind=GP)**6
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               r7(i,j,k) = (r2(i,j,k)*r6(i,j,k)-r5(i,j,k)*r3(i,j,k))*tmp
               r3(i,j,k) = (r3(i,j,k)*r4(i,j,k)-r6(i,j,k)*r1(i,j,k))*tmp
               r1(i,j,k) = (r1(i,j,k)*r5(i,j,k)-r4(i,j,k)*r2(i,j,k))*tmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,r7,d,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r3,e,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r1,f,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE aprodre3
