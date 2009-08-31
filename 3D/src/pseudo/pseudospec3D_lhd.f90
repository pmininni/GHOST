!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Extra subroutines to compute nonlinear terms in the 
! incompressible Leray-alpha model HD equations in 3D 
! using a pseudo-spectral method.  You should use the 
! FFTPLANS and MPIVARS modules (see the file 'fftp_mod.f90') 
! in each program that calls any of the subroutines in this 
! file.
!
! NOTATION: index 'i' is 'x'
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2007 Jonathan Pietarila Graham
!      National Center for Atmospheric Research.
!      e-mail: jgraham@ucar.edu
!=================================================================

!*****************************************************************
      SUBROUTINE aprodre3(a,b,c,d,e,f,alp)
!-----------------------------------------------------------------
!
! Three-dimensional product of the field As and 
! grad(A) in real space. The components of the 
! field A are given by the matrixes a,b and c
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     d  : product (As.grad)A_x in Fourier space [output]
!     e  : product (As.grad)A_y in Fourier space [output]
!     f  : product (As.grad)A_z in Fourier space [output]
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
      COMPLEX(KIND=GP), DIMENSION(n,n,ista:iend) :: c1,c2,c3
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)    :: r1,r2
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)    :: r3,r4
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend)    :: rx,ry,rz
      REAL(KIND=GP), INTENT(IN) :: alp
      REAL(KIND=GP)             :: tmp
      INTEGER          :: i,j,k

!
! Computes As
!
      CALL smooth3(a,b,c,c1,e,f,alp)
!
! Computes (As_x.dx)A_dir
!
      CALL derivk3(a,c2,1)
      CALL derivk3(b,c3,1)
      CALL derivk3(c,d,1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c3,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,d,r4,MPI_COMM_WORLD)

      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               rx(i,j,k) = r1(i,j,k)*r2(i,j,k)
               ry(i,j,k) = r1(i,j,k)*r3(i,j,k)
               rz(i,j,k) = r1(i,j,k)*r4(i,j,k)
            END DO
         END DO
      END DO
!
! Computes (As_y.dy)A_dir
!
      CALL derivk3(a,c2,2)
      CALL derivk3(b,c3,2)
      CALL derivk3(c,d,2)
      CALL fftp3d_complex_to_real(plancr,e,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c3,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,d,r4,MPI_COMM_WORLD)

      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               rx(i,j,k) = rx(i,j,k)+r1(i,j,k)*r2(i,j,k)
               ry(i,j,k) = ry(i,j,k)+r1(i,j,k)*r3(i,j,k)
               rz(i,j,k) = rz(i,j,k)+r1(i,j,k)*r4(i,j,k)
            END DO
         END DO
      END DO
!
! Computes (As_z.dz)A_dir
!
      CALL derivk3(a,c2,3)
      CALL derivk3(b,c3,3)
      CALL derivk3(c,d,3)
      CALL fftp3d_complex_to_real(plancr,f,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c3,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,d,r4,MPI_COMM_WORLD)

      tmp = 1./real(n,kind=GP)**6
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               rx(i,j,k) = (rx(i,j,k)+r1(i,j,k)*r2(i,j,k))*tmp
               ry(i,j,k) = (ry(i,j,k)+r1(i,j,k)*r3(i,j,k))*tmp
               rz(i,j,k) = (rz(i,j,k)+r1(i,j,k)*r4(i,j,k))*tmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,rx,d,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,ry,e,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,rz,f,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE aprodre3

!*****************************************************************
      SUBROUTINE aenergy(a,b,c,d,alp,kin)
!-----------------------------------------------------------------
!
! Computes the mean leray-alpha-energy of a vector 
! field. The output is valid only in the first 
! node.
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     d  : at the output contains the energy
!     alp: value of alpha
!     kin: =1 computes the kinetic energy
!          =0 computes the magnetic energy
!
      USE fprecision
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT) :: d
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      REAL(KIND=GP), INTENT(IN)    :: alp
      INTEGER, INTENT(IN) :: kin

      CALL energy(a,b,c,d,kin)

      RETURN
      END SUBROUTINE aenergy

!*****************************************************************
      SUBROUTINE across(a,b,c,d,e,f,g,alp,kin)
!-----------------------------------------------------------------
!
! Computes the leray-alpha inner product of two vector 
! fields. The output is valid only in the first 
! node.
!
! Parameters
!     a  : first field x-component
!     b  : first field y-component
!     c  : first field z-component
!     d  : second field x-component
!     e  : second field y-component
!     f  : second field z-component
!     g  : at the output contains the inner product
!     alp: value of alpha
!     kin: =1 computes the inner product of the fields
!          =0 computes the inner product of the curls
!
      USE fprecision
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT) :: g
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: d,e,f
      REAL(KIND=GP), INTENT(IN)    :: alp
      INTEGER, INTENT(IN) :: kin

      CALL cross(a,b,c,d,e,f,g,kin)

      RETURN
      END SUBROUTINE across
