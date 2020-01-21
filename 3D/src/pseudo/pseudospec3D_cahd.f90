!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Extra subroutines to compute nonlinear terms in the 
! incompressible Clark-alpha HD equations in 3D using a 
! pseudo-spectral method. You should use the FFTPLANS 
! and MPIVARS modules (see the file 'fftp_mod.f90') in 
! each program that call any of the subroutines in this 
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
! (1-.5 alpha^2 laplacian)w_s x u -.5 alpha^2(w_s x (laplacian u)
! + (laplacian w_s) x u ) carried out in real space. The Fourier 
! components of the field v are given by the matrixes a,b and c. 
! Fourier output in d, e, and f.
!
! Parameters
!     a  : input matrix in the x-direction in Fourier space
!     b  : input matrix in the y-direction in Fourier space
!     c  : input matrix in the z-direction in Fourier space
!     d  : output _x in Fourier space
!     e  : output _y in Fourier space
!     f  : output _z in Fourier space
!     alp: value of alpha
!
      USE fprecision
      USE kes
      USE commtypes
      USE mpivars
      USE grid
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2,c3
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c4,c5,c6
      REAL(KIND=GP),    DIMENSION(nx,ny,ksta:kend) :: r1,r2,r3
      REAL(KIND=GP),    DIMENSION(nx,ny,ksta:kend) :: r4,r5,r6
      REAL(KIND=GP),    DIMENSION(nx,ny,ksta:kend) :: r7,r8,r9
      REAL(KIND=GP),    INTENT(IN) :: alp
      REAL(KIND=GP)                :: tmp,rx,ry
      INTEGER                      :: i,j,k

! v -> u: c_{1-3}
      CALL smooth3(a,b,c,c1,c2,c3,alp)

! u -> w_s: c_{4-6}
      CALL rotor3(c2,c3,c4,1)
      CALL rotor3(c1,c3,c5,2)
      CALL rotor3(c1,c2,c6,3)

! r_{1-3} = w_s * n^3
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               d(k,j,i) = c4(k,j,i)
               e(k,j,i) = c5(k,j,i)
               f(k,j,i) = c6(k,j,i)
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,d,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,e,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,f,r3,MPI_COMM_WORLD)

! r_{4-6} = u * n^3
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               d(k,j,i) = c1(k,j,i)
               e(k,j,i) = c2(k,j,i)
               f(k,j,i) = c3(k,j,i)
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,d,r4,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,e,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,f,r6,MPI_COMM_WORLD)

! r_{7-9} = w_s x u
      tmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r7(i,j,k) = (r2(i,j,k)*r6(i,j,k)-r5(i,j,k)*r3(i,j,k))*tmp
               r8(i,j,k) = (r3(i,j,k)*r4(i,j,k)-r6(i,j,k)*r1(i,j,k))*tmp
               r9(i,j,k) = (r1(i,j,k)*r5(i,j,k)-r4(i,j,k)*r2(i,j,k))*tmp
            END DO
         END DO
      END DO

! d-f = w_s x u
      CALL fftp3d_real_to_complex(planrc,r7,d,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r8,e,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r9,f,MPI_COMM_WORLD)

! d-f = w_s x u - .5 alpha^2 \nabla^2(w_s x u)
      tmp = 0.5 * alp**2
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               d(k,j,i) = (1+kk2(k,j,i)*tmp)*d(k,j,i)
            END DO
         END DO
      END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               e(k,j,i) = (1+kk2(k,j,i)*tmp)*e(k,j,i)
            END DO
         END DO
      END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               f(k,j,i) = (1+kk2(k,j,i)*tmp)*f(k,j,i)
            END DO
         END DO
      END DO

! c_{4-6} = \nabla^2 w_s
      CALL laplak3(c4,c4)
      CALL laplak3(c5,c5)
      CALL laplak3(c6,c6)

! r_{7-9} = \nalba^2 w_s * n^3
! c_{4-6} DESTROYED
      CALL fftp3d_complex_to_real(plancr,c4,r7,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c5,r8,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c6,r9,MPI_COMM_WORLD)

! r_{7-9} = - .5 \alpha^2 \nabla^2 w_s x u
      tmp = (-0.5 * alp**2)/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               rx = (r8(i,j,k)*r6(i,j,k)-r5(i,j,k)*r9(i,j,k))*tmp
               ry = (r9(i,j,k)*r4(i,j,k)-r6(i,j,k)*r7(i,j,k))*tmp
               r9(i,j,k) = (r7(i,j,k)*r5(i,j,k)-r4(i,j,k)*r8(i,j,k))*tmp
               r7(i,j,k) = rx
               r8(i,j,k) = ry
            END DO
         END DO
      END DO

! c_{1-3} = \nabla^2 u
      CALL laplak3(c1,c1)
      CALL laplak3(c2,c2)
      CALL laplak3(c3,c3)

! r_{4-6} = \nalba^2 u * n^3
! c_{1-3} DESTROYED
      CALL fftp3d_complex_to_real(plancr,c1,r4,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c3,r6,MPI_COMM_WORLD)

! r_{x-z} = - .5 \alpha^2 \nabla^2 w_s x u- .5 \alpha^2 w_s x \nabla^2 u
      tmp = (0.5 * alp**2)/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r7(i,j,k)=r7(i,j,k)-(r2(i,j,k)*r6(i,j,k)-r5(i,j,k) &
                         *r3(i,j,k))*tmp
               r8(i,j,k)=r8(i,j,k)-(r3(i,j,k)*r4(i,j,k)-r6(i,j,k) &
                         *r1(i,j,k))*tmp
               r9(i,j,k)=r9(i,j,k)-(r1(i,j,k)*r5(i,j,k)-r4(i,j,k) &
                         *r2(i,j,k))*tmp
            END DO
         END DO
      END DO

! c_{1-3} = - .5 \alpha^2 \nabla^2 w_s x u- .5 \alpha^2 w_s x \nabla^2 u
      CALL fftp3d_real_to_complex(planrc,r7,c1,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r8,c2,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r9,c3,MPI_COMM_WORLD)

! d-f = w_s x u - .5 \alpha^2 \nabla^2(w_s x u)
! - .5 \alpha^2 \nabla^2 w_s x u- .5 \alpha^2 w_s x \nabla^2 u
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               d(k,j,i) = d(k,j,i) + c1(k,j,i)
            END DO
         END DO
      END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               e(k,j,i) = e(k,j,i) + c2(k,j,i)
            END DO
         END DO
      END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               f(k,j,i) = f(k,j,i) + c3(k,j,i)
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE aprodre3
