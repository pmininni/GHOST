!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines to compute spatial derivatives and nonlinear 
! terms in Navier-Stokes, MHD and Hall-MHD equations in 3D 
! using a pseudo-spectral method. You should use the FFTPLANS 
! and MPIVARS modules (see the file 'fftp_mod.f90') in each 
! program that calls any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!      
! 2003 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!=================================================================

!*****************************************************************
      SUBROUTINE derivk3(a,b,dir)
!-----------------------------------------------------------------
!
! Three-dimensional derivative of the matrix 'a'
!
! Parameters
!     a  : input matrix
!     b  : at the output contains the derivative da/dk_dir
!     dir: =1 derivative in the x-direction
!          =2 derivative in the y-direction
!          =3 derivative in the z-direction
!
      USE kes
      USE var
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: b
      INTEGER, INTENT(IN) :: dir
      INTEGER             :: i,j,k

!
! Derivative in the x-direction
!
      IF (dir.eq.1) THEN
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  b(k,j,i) = im*kx(i)*a(k,j,i)
               END DO
            END DO
         END DO
!
! Derivative in the y-direction
!
      ELSE IF (dir.eq.2) THEN
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  b(k,j,i) = im*ky(j)*a(k,j,i)
               END DO
            END DO
         END DO
!
! Derivative in the z-direction
!
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  b(k,j,i) = im*kz(k)*a(k,j,i)
               END DO
            END DO
         END DO
      ENDIF

      RETURN
      END SUBROUTINE derivk3

!*****************************************************************
      SUBROUTINE laplak3(a,b)
!-----------------------------------------------------------------
!
! Three-dimensional Laplacian of the matrix 'a'
!
! Parameters
!     a: input matrix
!     b: at the output contains the Laplacian
!
      USE kes
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: b
      INTEGER :: i,j,k

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               b(k,j,i) = -kk2(k,j,i)*a(k,j,i)
            END DO
         END DO
      END DO
      
      RETURN
      END SUBROUTINE laplak3

!*****************************************************************
      SUBROUTINE rotor3(a,b,c,dir)
!-----------------------------------------------------------------
!
! Computes the curl of the vector field A in Fourier
! space. The needed components of the field A are 
! given by the matrixes a and b, and the order must 
! follow the right hand convention.
!
! Parameters
!     a  : input matrix
!     b  : input matrix
!     c  : at the output contains curl(A)_dir
!     dir: =1 computes the x-component
!          =2 computes the y-component
!          =3 computes the z-component
!
      USE fprecision
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: c
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2
      INTEGER, INTENT(IN) :: dir
      INTEGER             :: i,j,k

!
! Computes the x-component
!
      IF (dir.eq.1) THEN
         CALL derivk3(a,c1,3)
         CALL derivk3(b,c2,2)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  c(k,j,i) = c2(k,j,i)-c1(k,j,i)
               END DO
            END DO
         END DO
!
! Computes the y-component
!
      ELSE IF (dir.eq.2) THEN
         CALL derivk3(a,c1,3)
         CALL derivk3(b,c2,1)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  c(k,j,i) = c1(k,j,i)-c2(k,j,i)
               END DO
            END DO
         END DO
!
! Computes the z-component
!
      ELSE
         CALL derivk3(a,c1,2)
         CALL derivk3(b,c2,1)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  c(k,j,i) = c2(k,j,i)-c1(k,j,i)
               END DO
            END DO
         END DO
      ENDIF

      RETURN
      END SUBROUTINE rotor3

!*****************************************************************
      SUBROUTINE gradre3(a,b,c,d,e,f)
!-----------------------------------------------------------------
!
! Three-dimensional inner product A.grad(A) in 
! real space. The components of the field A are 
! given by the matrixes a, b and c
!
! Parameters
!     a: input matrix in the x-direction
!     b: input matrix in the y-direction
!     c: input matrix in the z-direction
!     d: product (A.grad)A_x in Fourier space [output]
!     e: product (A.grad)A_y in Fourier space [output]
!     f: product (A.grad)A_z in Fourier space [output]
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c3,c4
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r1,r2
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r3,r4
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: rx,ry,rz
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k

!
! Computes (A_x.dx)A_dir
!
      c1 = a
      CALL derivk3(a,c2,1)
      CALL derivk3(b,c3,1)
      CALL derivk3(c,c4,1)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c3,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c4,r4,MPI_COMM_WORLD)

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               rx(i,j,k) = r1(i,j,k)*r2(i,j,k)
               ry(i,j,k) = r1(i,j,k)*r3(i,j,k)
               rz(i,j,k) = r1(i,j,k)*r4(i,j,k)
            END DO
         END DO
      END DO
!
! Computes (A_y.dy)A_dir
!
      c1 = b
      CALL derivk3(a,c2,2)
      CALL derivk3(b,c3,2)
      CALL derivk3(c,c4,2)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c3,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c4,r4,MPI_COMM_WORLD)

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               rx(i,j,k) = rx(i,j,k)+r1(i,j,k)*r2(i,j,k)
               ry(i,j,k) = ry(i,j,k)+r1(i,j,k)*r3(i,j,k)
               rz(i,j,k) = rz(i,j,k)+r1(i,j,k)*r4(i,j,k)
            END DO
         END DO
      END DO
!
! Computes (A_z.dz)A_dir
!
      c1 = c
      CALL derivk3(a,c2,3)
      CALL derivk3(b,c3,3)
      CALL derivk3(c,c4,3)
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c3,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c4,r4,MPI_COMM_WORLD)

      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
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
      END SUBROUTINE gradre3

!*****************************************************************
      SUBROUTINE prodre3(a,b,c,d,e,f)
!-----------------------------------------------------------------
!
! Three-dimensional cross product curl(A)xA in 
! real space. The components of the field A are 
! given by the matrixes a, b and c
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     d  : product [curl(A)xA]_x in Fourier space [output]
!     e  : product [curl(A)xA]_y in Fourier space [output]
!     f  : product [curl(A)xA]_z in Fourier space [output]
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: d,e,f
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r1,r2
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r3,r4
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r5,r6
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)    :: r7
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k

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
! Computes A
!
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               d(k,j,i) = a(k,j,i)
               e(k,j,i) = b(k,j,i)
               f(k,j,i) = c(k,j,i)
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,d,r4,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,e,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,f,r6,MPI_COMM_WORLD)
!
! Computes curl(A)xA
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
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
      END SUBROUTINE prodre3

!*****************************************************************
      SUBROUTINE nonlhd3(a,b,c,g,dir)
!-----------------------------------------------------------------
!
! Computes the nonlinear terms in 3D Navier-Stokes 
! equation. It takes the components of (v.grad)v or 
! curl(v)xv as input matrixes and computes 
! -(v.grad)v-grad(p) or -curl(v)xv-grad(p) in 
! Fourier space, with the pressure chosen to satisfy 
! the incompressibility condition.
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     g  : at the output contains the result in Fourier space
!     dir: =1 computes the x-component
!          =2 computes the y-component
!          =3 computes the z-component
!
      USE fprecision
      USE kes
      USE grid
      USE mpivars
      USE commtypes
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: g
      INTEGER, INTENT(IN) :: dir
      INTEGER             :: i,j,k

!
! Computes the x-component
!
      IF (dir.eq.1) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k)
            DO j = 1,ny
               DO k = 1,nz
                  g(k,j,1) = -a(k,j,1)
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k)
               DO j = 1,ny
                  DO k = 1,nz
                     g(k,j,i) = -a(k,j,i)+kx(i)*(kx(i)*a(k,j,i) &
                       +ky(j)*b(k,j,i)+kz(k)*c(k,j,i))/kk2(k,j,i)
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
               DO j = 1,ny
                  DO k = 1,nz
                     g(k,j,i) = -a(k,j,i)+kx(i)*(kx(i)*a(k,j,i) &
                       +ky(j)*b(k,j,i)+kz(k)*c(k,j,i))/kk2(k,j,i)
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the y-component
!
      ELSE IF (dir.eq.2) THEN
!$omp parallel do if (iend-ista.ge.nth) private (k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth)
            DO k = 1,nz
               g(k,1,i) = -b(k,1,i)
            END DO
         END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 2,ny
               DO k = 1,nz
                  g(k,j,i) = -b(k,j,i)+ky(j)*(kx(i)*a(k,j,i) &
                    +ky(j)*b(k,j,i)+kz(k)*c(k,j,i))/kk2(k,j,i)
               END DO
            END DO
         END DO
!
! Computes the z-component
!
      ELSE

!$omp parallel do if (iend-ista.ge.nth) private (j)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth)
            DO j = 1,ny
               g(1,j,i) = -c(1,j,i)
            END DO
         END DO
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 2,nz
                  g(k,j,i) = -c(k,j,i)+kz(k)*(kx(i)*a(k,j,i) &
                    +ky(j)*b(k,j,i)+kz(k)*c(k,j,i))/kk2(k,j,i)
               END DO
            END DO
         END DO
      ENDIF
      RETURN
      END SUBROUTINE nonlhd3

!*****************************************************************
      SUBROUTINE energy(a,b,c,d,kin)
!-----------------------------------------------------------------
!
! Computes the mean energy of a vector field.
! The output is only valid in the first node.
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     d  : at the output contains the energy
!     kin: =0 computes the magnetic energy
!          =1 computes the kinetic energy
!          =2 computes the magnetic enstrophy
!
      USE fprecision
      USE commtypes
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)          :: c1,c2,c3
      DOUBLE PRECISION, INTENT(OUT) :: d
      DOUBLE PRECISION              :: dloc
      REAL(KIND=GP)                 :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k

      dloc = 0.0D0
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2

!
! Computes the kinetic energy
!
      IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz
                  dloc = dloc+(abs(a(k,j,1))**2+abs(b(k,j,1))**2+ &
                         abs(c(k,j,1))**2)*tmp
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:dloc)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:dloc)
               DO j = 1,ny
                  DO k = 1,nz
                     dloc = dloc+2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+ &
                            abs(c(k,j,i))**2)*tmp
                  END DO
               END DO
            END DO
          ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:dloc)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:dloc)
               DO j = 1,ny
                  DO k = 1,nz
                     dloc = dloc+2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+ &
                            abs(c(k,j,i))**2)*tmp
                  END DO
               END DO
            END DO
          ENDIF
!
! Computes the magnetic energy
!
      ELSE IF (kin.eq.0) THEN
         CALL rotor3(b,c,c1,1)
         CALL rotor3(a,c,c2,2)
         CALL rotor3(a,b,c3,3)
         IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz
                  dloc = dloc+(abs(c1(k,j,1))**2+abs(c2(k,j,1))**2+ &
                         abs(c3(k,j,1))**2)*tmp
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:dloc)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:dloc)
               DO j = 1,ny
                  DO k = 1,nz
                     dloc = dloc+2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+ &
                            abs(c3(k,j,i))**2)*tmp
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:dloc)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:dloc)
               DO j = 1,ny
                  DO k = 1,nz
                     dloc = dloc+2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+ &
                            abs(c3(k,j,i))**2)*tmp
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the magnetic enstrophy
!
      ELSE
         CALL laplak3(a,c1)
         CALL laplak3(b,c2)
         CALL laplak3(c,c3)
         IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz
                  dloc = dloc+(abs(c1(k,j,1))**2+abs(c2(k,j,1))**2+ &
                         abs(c3(k,j,1))**2)*tmp
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:dloc)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:dloc)
               DO j = 1,ny
                  DO k = 1,nz
                     dloc = dloc+2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+ &
                            abs(c3(k,j,i))**2)*tmp
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:dloc)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:dloc)
               DO j = 1,ny
                  DO k = 1,nz
                     dloc = dloc+2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+ &
                            abs(c3(k,j,i))**2)*tmp
                  END DO
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(dloc,d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE energy

!*****************************************************************
      SUBROUTINE helicity(a,b,c,d)
!-----------------------------------------------------------------
!
! Computes the mean helicity of a vector field.
! The output is only valid in the first node.
!
! Parameters
!     a: input matrix in the x-direction
!     b: input matrix in the y-direction
!     c: input matrix in the z-direction
!     d: at the output contains the helicity
!
      USE fprecision
      USE commtypes
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)             :: c1
      DOUBLE PRECISION, INTENT(OUT) :: d
      DOUBLE PRECISION              :: dloc
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k

      dloc = 0.0D0
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2

      CALL rotor3(b,c,c1,1)
      IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:dloc)
         DO j = 1,ny
            DO k = 1,nz
                 dloc = dloc+real(a(k,j,1)*conjg(c1(k,j,1)))*tmp
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:dloc)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz
                  dloc = dloc+2*real(a(k,j,i)*conjg(c1(k,j,i)))*tmp
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:dloc)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz
                  dloc = dloc+2*real(a(k,j,i)*conjg(c1(k,j,i)))*tmp
               END DO
            END DO
         END DO
      ENDIF
      CALL rotor3(a,c,c1,2)
      IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:dloc)
         DO j = 1,ny
            DO k = 1,nz
               dloc = dloc+real(b(k,j,1)*conjg(c1(k,j,1)))*tmp
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:dloc)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz
                  dloc = dloc+2*real(b(k,j,i)*conjg(c1(k,j,i)))*tmp
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:dloc)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz
                  dloc = dloc+2*real(b(k,j,i)*conjg(c1(k,j,i)))*tmp
               END DO
            END DO
         END DO
      ENDIF
      CALL rotor3(a,b,c1,3)
      IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:dloc)
         DO j = 1,ny
            DO k = 1,nz
               dloc = dloc+real(c(k,j,1)*conjg(c1(k,j,1)))*tmp
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:dloc)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz
                  dloc = dloc+2*real(c(k,j,i)*conjg(c1(k,j,i)))*tmp
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:dloc)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz
                  dloc = dloc+2*real(c(k,j,i)*conjg(c1(k,j,i)))*tmp
               END DO
            END DO
         END DO
      ENDIF

      CALL MPI_REDUCE(dloc,d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE helicity

!*****************************************************************
      SUBROUTINE cross(a,b,c,d,e,f,g,kin)
!-----------------------------------------------------------------
!
! Computes the cross helicity or averaged inner 
! product of two vector fields. The output is 
! only valid in the first node.
!
! Parameters
!     a  : first field x-component
!     b  : first field y-component
!     c  : first field z-component
!     d  : second field x-component
!     e  : second field y-component
!     f  : second field z-component
!     g  : at the output contains the inner product
!     kin: =0 computes the inner product of the curls
!          =1 computes the inner product of the fields
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: d,e,f
      DOUBLE PRECISION, INTENT(OUT) :: g
      DOUBLE PRECISION              :: gloc
      REAL(KIND=GP)                 :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k

      gloc = 0.0D0
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!
! Computes the averaged inner product between the fields
!
      IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:gloc)
            DO j = 1,ny
               DO k = 1,nz
                  gloc = gloc+real(a(k,j,1)*conjg(d(k,j,1))+      &
                        b(k,j,1)*conjg(e(k,j,1))+c(k,j,1)*        &
                        conjg(f(k,j,1)))*tmp
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:gloc)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:gloc)
               DO j = 1,ny
                  DO k = 1,nz
                     gloc = gloc+2*real(a(k,j,i)*conjg(d(k,j,i))+ &
                           b(k,j,i)*conjg(e(k,j,i))+c(k,j,i)*     &
                           conjg(f(k,j,i)))*tmp
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:gloc)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:gloc)
               DO j = 1,ny
                  DO k = 1,nz
                     gloc = gloc+2*real(a(k,j,i)*conjg(d(k,j,i))+ &
                           b(k,j,i)*conjg(e(k,j,i))+c(k,j,i)*     &
                           conjg(f(k,j,i)))*tmp
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the averaged inner product between 
! the curl of the fields
!
      ELSE
         IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:gloc)
            DO j = 1,ny
               DO k = 1,nz
                  gloc = gloc+real(a(k,j,1)*conjg(d(k,j,1))+      &
                        b(k,j,1)*conjg(e(k,j,1))+c(k,j,1)*        &
                        conjg(f(k,j,1)))*kk2(k,j,1)*tmp
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:gloc)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:gloc)
               DO j = 1,ny
                  DO k = 1,nz
                     gloc = gloc+2*real(a(k,j,i)*conjg(d(k,j,i))+ &
                           b(k,j,i)*conjg(e(k,j,i))+c(k,j,i)*     &
                           conjg(f(k,j,i)))*kk2(k,j,i)*tmp
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:gloc)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:gloc)
               DO j = 1,ny
                  DO k = 1,nz
                     gloc = gloc+2*real(a(k,j,i)*conjg(d(k,j,i))+ &
                           b(k,j,i)*conjg(e(k,j,i))+c(k,j,i)*     &
                           conjg(f(k,j,i)))*kk2(k,j,i)*tmp
                  END DO
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(gloc,g,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE cross

!*****************************************************************
      SUBROUTINE maxabs(a,b,c,d,kin)
!-----------------------------------------------------------------
!
! Computes the maximum absolute value of the field 
! vorticity, current density, or of the original field. 
! The output is only valid in the first node.
!
! Parameters
!     a  : field x-component
!     b  : field y-component
!     c  : field z-component
!     d  : at the output contains the maximum value
!     kin: =0 computes the maximum of vorticity
!          =1 computes the maximum of current density
!          =2 computes the maximum of the field
!
      USE fprecision
      USE commtypes
      USE fft
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)          :: c1,c2,c3
      REAL(KIND=GP), DIMENSION(nx,ny,ksta:kend)             :: r1,r2,r3
      REAL(KIND=GP), INTENT(OUT)   :: d
      REAL(KIND=GP)                :: dloc
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k

      IF (kin.eq.0) THEN
         CALL rotor3(b,c,c1,1)
         CALL rotor3(a,c,c2,2)
         CALL rotor3(a,b,c3,3)
      ELSE IF (kin.eq.1) THEN
         CALL laplak3(a,c1)
         CALL laplak3(b,c2)
         CALL laplak3(c,c3)
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  c1(k,j,i) = a(k,j,i)
                  c2(k,j,i) = b(k,j,i)
                  c3(k,j,i) = c(k,j,i)
               END DO
            END DO
         END DO
      ENDIF
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c3,r3,MPI_COMM_WORLD)
      dloc = 0.0_GP
!$omp parallel do if (kend-ksta.ge.nth) private (j,i) reduction(max:dloc)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i) reduction(max:dloc)
         DO j = 1,ny
            DO i = 1,nx
               dloc = max(dloc, &
                          sqrt(r1(i,j,k)**2+r2(i,j,k)**2+r3(i,j,k)**2))
            END DO
         END DO
      END DO
      dloc = dloc/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
      CALL MPI_REDUCE(dloc,d,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE maxabs

!*****************************************************************
      SUBROUTINE hdcheck(a,b,c,d,e,f,t,dt,hel,chk)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of energy, 
! helicity, and null divergency of the velocity field
!
! Output files contain:
! 'balance.txt':  time, <v^2>, <omega^2>, mechanic injection rate
! 'helicity.txt': time, kinetic helicity
! 'divergence.txt' [OPTIONAL]: time, <(div.v)^2>
!
! Parameters
!     a  : velocity field in the x-direction
!     b  : velocity field in the y-direction
!     c  : velocity field in the z-direction
!     d  : force in the x-direction
!     e  : force in the y-direction
!     f  : force in the z-direction
!     t  : number of time steps made
!     dt : time step
!     hel: =0 skips kinetic helicity computation
!          =1 computes the kinetic helicity
!     chk: =0 skips divergency check
!          =1 performs divergency check
!
      USE fprecision
      USE commtypes
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)          :: c1,c2,c3
      DOUBLE PRECISION    :: eng,ens,pot,khe
      DOUBLE PRECISION    :: div,tmp
      REAL(KIND=GP)       :: dt
      REAL(KIND=GP)       :: tmq
      INTEGER, INTENT(IN) :: hel,chk
      INTEGER, INTENT(IN) :: t
      INTEGER             :: i,j,k

      div = 0.0D0
      tmp = 0.0D0
      tmq = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2

!
! Computes the mean square value of
! the divergence of the vector field
!
      IF (chk.eq.1) THEN

      CALL derivk3(a,c1,1)
      CALL derivk3(b,c2,2)
      CALL derivk3(c,c3,3)
      IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:tmp)
         DO j = 1,ny
            DO k = 1,nz
               tmp = tmp+abs(c1(k,j,1)+c2(k,j,1)+c3(k,j,1))**2*tmq
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:tmp)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:tmp)
            DO j = 1,ny
               DO k = 1,nz
                  tmp = tmp+2*abs(c1(k,j,i)+c2(k,j,i)+c3(k,j,i))**2*tmq
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:tmp)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:tmp)
            DO j = 1,ny
               DO k = 1,nz
                  tmp = tmp+2*abs(c1(k,j,i)+c2(k,j,i)+c3(k,j,i))**2*tmq
               END DO
            END DO
         END DO
      ENDIF
      CALL MPI_REDUCE(tmp,div,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      ENDIF
!
! Computes the mean energy, enstrophy, and kinetic helicity
!
      CALL energy(a,b,c,eng,1)
      CALL energy(a,b,c,ens,0)
      IF (hel.eq.1) THEN
         CALL helicity(a,b,c,khe)
      ENDIF
!
! Computes the energy injection rate
!
      CALL cross(a,b,c,d,e,f,pot,1)
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='balance.txt',position='append')
         WRITE(1,10) (t-1)*dt,eng,ens,pot
   10    FORMAT( E13.6,E25.18,E25.18,E25.18 )
         CLOSE(1)
         IF (hel.eq.1) THEN
            OPEN(1,file='helicity.txt',position='append')
            WRITE(1,*) (t-1)*dt,khe
            CLOSE(1)
         ENDIF
         IF (chk.eq.1) THEN
            OPEN(1,file='divergence.txt',position='append')
            WRITE(1,*) (t-1)*dt,div
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE hdcheck

!*****************************************************************
      SUBROUTINE spectrum(a,b,c,nmb,kin,hel)
!-----------------------------------------------------------------
!
! Computes the energy and helicity power 
! spectrum. The output is written to a 
! file by the first node.
!
! Output files contain:
! 'kspectrum.XXX.txt': k, Ev(k)
! 'mspectrum.XXX.txt': k, Eb(k)
! 'khelicity.XXX.txt': k, Hv(k) (kinetic helicity spectrum)
! 'mhelicity.XXX.txt': k, Hb(k) (magnetic helicity spectrum)
! 'mhelicity.XXX.txt': k, Hb(k) (magnetic helicity spectrum)
! 'ghelicity.XXX.txt': k, G(k)  (generalized helicity in Hall-MHD)
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     nmb: the extension used when writting the file
!     kin: =0 computes the magnetic spectrum
!          =1 computes the kinetic spectrum
!          =2 skips energy spectrum computation
!     hel: =0 skips helicity spectrum computation
!          =1 computes the helicity spectrum
!
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nmax/2+1) :: Ek,Hk
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      INTEGER, INTENT(IN)          :: kin,hel
      INTEGER                      :: i
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Computes the energy and/or helicity spectra
!
      CALL spectrumc(a,b,c,kin,hel,Ek,Hk)
!
! Exports the energy spectrum to a file
!
      IF (kin.le.1) THEN
         IF (myrank.eq.0) THEN
            IF (kin.eq.1) THEN
               OPEN(1,file='kspectrum.' // nmb // '.txt')
            ELSE
               OPEN(1,file='mspectrum.' // nmb // '.txt')
            ENDIF
            DO i=1,nmax/2+1
               WRITE(1,FMT='(E13.6,E23.15)')  Dkk*i,Ek(i)
            END DO
            CLOSE(1)
         ENDIF
      ENDIF
!
! Exports the helicity spectrum to a file
!
      IF (hel.eq.1) THEN
         IF (myrank.eq.0) THEN
            IF (kin.eq.1) THEN
               OPEN(1,file='khelicity.' // nmb // '.txt')
            ELSE IF (kin.eq.0) THEN
               OPEN(1,file='mhelicity.' // nmb // '.txt')
            ELSE
               OPEN(1,file='ghelicity.' // nmb // '.txt')
            ENDIF
            DO i=1,nmax/2+1
               WRITE(1,FMT='(E13.6,E23.15)')  Dkk*i,Hk(i)
            END DO
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE spectrum

!*****************************************************************
      SUBROUTINE spectrumc(a,b,c,kin,hel,Ektot,Hktot)
!-----------------------------------------------------------------
!
! Computes the energy and helicity power 
! spectra, returning them.
!
! Parameters
!     a    : input matrix in the x-direction
!     b    : input matrix in the y-direction
!     c    : input matrix in the z-direction
!     kin  : =0 computes the magnetic spectrum
!            =1 computes the kinetic spectrum
!            =2 skips energy spectrum computation
!     hel  : =0 skips helicity spectrum computation
!            =1 computes the helicity spectrum
!     Ektot: output energy spectrum
!     Hktot: output helicity spectrum
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nmax/2+1)              :: Ek
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(nmax/2+1) :: Ektot, Hktot
      DOUBLE PRECISION    :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)          :: c1,c2,c3
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: kin,hel
      INTEGER             :: i,j,k
      INTEGER             :: kmn

!
! Computes the curl of the field if needed
!
      IF ((kin.eq.0).or.(hel.eq.1)) THEN
         CALL rotor3(b,c,c1,1)
         CALL rotor3(a,c,c2,2)
         CALL rotor3(a,b,c3,3)
      ENDIF
!
! Computes the kinetic energy spectrum
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      IF (kin.eq.1) THEN
         DO i = 1,nmax/2+1
            Ek   (i) = 0.0D0
            Ektot(i) = 0.0D0
         END DO
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(sqrt(kk2(k,j,1))/Dkk+.501)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = (abs(a(k,j,1))**2+abs(b(k,j,1))**2+        &
                            abs(c(k,j,1))**2)*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(sqrt(kk2(k,j,i))/Dkk+.501)
                     IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                        tmq = 2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+   &
                                 abs(c(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                  END DO
               END DO
            END DO
          ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(sqrt(kk2(k,j,i))/Dkk+.501)
                     IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                        tmq = 2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+   &
                                 abs(c(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                  END DO
               END DO
            END DO
          ENDIF
!
! Computes the magnetic energy spectrum
!
      ELSE IF (kin.eq.0) THEN
         DO i = 1,nmax/2+1
            Ek   (i) = 0.0D0
            Ektot(i) = 0.0D0
         END DO
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(sqrt(kk2(k,j,1))/Dkk+.501)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = (abs(c1(k,j,1))**2+abs(c2(k,j,1))**2+      &
                            abs(c3(k,j,1))**2)*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(sqrt(kk2(k,j,i))/Dkk+.501)
                     IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                        tmq = 2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+ &
                                 abs(c3(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(sqrt(kk2(k,j,i))/Dkk+.501)
                     IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                        tmq = 2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+ &
                                 abs(c3(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
!
      IF (kin.le.1) THEN
         CALL MPI_ALLREDUCE(Ek,Ektot,nmax/2+1,MPI_DOUBLE_PRECISION,   &
                         MPI_SUM,MPI_COMM_WORLD,ierr)
      END IF
!
! Computes the helicity spectrum
!
      IF (hel.eq.1) THEN
         DO i = 1,nmax/2+1
            Ek(i) = 0.0D0
            Hktot(i) = 0.0D0
         END DO
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(sqrt(kk2(k,j,1))/Dkk+.501)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = (real(a(k,j,1)*conjg(c1(k,j,1)))+          &
                            real(b(k,j,1)*conjg(c2(k,j,1)))+          &
                            real(c(k,j,1)*conjg(c3(k,j,1))))*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(sqrt(kk2(k,j,i))/Dkk+.501)
                     IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+     &
                                 real(b(k,j,i)*conjg(c2(k,j,i)))+     &
                                 real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                 END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(sqrt(kk2(k,j,i))/Dkk+.501)
                     IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+     &
                                 real(b(k,j,i)*conjg(c2(k,j,i)))+     &
                                 real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the reduction between nodes
!
         CALL MPI_ALLREDUCE(Ek,Hktot,nmax/2+1,MPI_DOUBLE_PRECISION,   &
                         MPI_SUM,MPI_COMM_WORLD,ierr)
      ENDIF

      RETURN
      END SUBROUTINE spectrumc

!*****************************************************************
      SUBROUTINE spectr1d(a,nmb,cmp,dir)
!-----------------------------------------------------------------
!
! Computes the 1D longitudinal or transverse kinetic energy 
! spectrum following the nomenclature of Monin and Yaglom. 
! The k-shells are planes with normal (0,0,k_dir), with 
! k_dir=(0,...,nmax/2)*dk. The output is written to a file
! by the first node.
!
! Output files contain:
! 'kspec1dIJ.XXX.txt' [I=x,y,z (field component)][J=x,y,z (dir)]:
!     k_dir, E_comp(k_dir)
!
! Parameters
!     a  : input matrix with a field component
!     nmb: the extension used when writting the file
!     cmp: =1 the input matrix is v_x
!          =2 the input matrix is v_y
!          =3 the input matrix is v_z
!     dir: =1 computes the 1D spectrum in k_x
!          =2 computes the 1D spectrum in k_y
!          =3 computes the 1D spectrum in k_z
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nmax/2+1) :: Ek,Ektot
      DOUBLE PRECISION    :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a
      REAL(KIND=GP)       :: Dn(3)
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: cmp,dir
      INTEGER             :: i,j,k
      INTEGER             :: n(3)
      INTEGER             :: kmn
      CHARACTER(len=3)    :: coord
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,nmax/2+1
         Ek(i) = 0.
      END DO
!
! Computes the kinetic energy spectrum
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      IF (dir.eq.1) THEN ! E(k_x)
         IF (ista.eq.1) THEN
            kmn = int(abs(kx(1))*Lx+1)
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  tmq = abs(a(k,j,1))**2*tmp
!$omp atomic
                  Ek(kmn) = Ek(kmn)+tmq                       
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
               kmn = int(abs(kx(i))*Lx+1)
               IF (kmn.le.nx/2+1) THEN
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
                  DO j = 1,ny
                     DO k = 1,nz
                        tmq = 2*abs(a(k,j,i))**2*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  END DO
               ENDIF
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
               kmn = int(abs(kx(i))*Lx+1)
               IF (kmn.le.nx/2+1) THEN
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
                  DO j = 1,ny
                     DO k = 1,nz
                        tmq = 2*abs(a(k,j,i))**2*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  END DO
               ENDIF
            END DO
         ENDIF
      ELSEIF (dir.eq.2) THEN ! E(k_y)
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,ny
               kmn = int(abs(ky(j))*Ly+1)
               IF (kmn.le.ny/2+1) THEN
                  DO k = 1,nz
                     tmq = abs(a(k,j,1))**2*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq                       
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  kmn = int(abs(ky(j))*Ly+1)
                  IF (kmn.le.ny/2+1) THEN
                     DO k = 1,nz
                        tmq = 2*abs(a(k,j,i))**2*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  kmn = int(abs(ky(j))*Ly+1)
                  IF (kmn.le.ny/2+1) THEN
                     DO k = 1,nz
                        tmq = 2*abs(a(k,j,i))**2*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  ENDIF
               END DO
            END DO
         ENDIF
      ELSE                   ! E(k_z)
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(abs(kz(k))*Lz+1)
                  IF (kmn.le.nz/2+1) THEN
                     tmq = abs(a(k,j,1))**2*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq                       
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(abs(kz(k))*Lz+1)
                     IF (kmn.le.nz/2+1) THEN
                        tmq = 2*abs(a(k,j,i))**2*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(abs(kz(k))*Lz+1)
                     IF (kmn.le.nz/2+1) THEN
                        tmq = 2*abs(a(k,j,i))**2*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Ek,Ektot,nmax/2+1,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         coord = 'xyz'
         n  = (/nx,ny,nz/)
         Dn = (/Dkx,Dky,Dkz/)
         OPEN(1,file='kspec1d' // coord(cmp:cmp) // coord(dir:dir)  &
              // '.' // nmb // '.txt')
         DO i=1,n(dir)/2+1
            WRITE(1,FMT='(E13.6,E23.15)')  Dn(dir)*(i-1),Ektot(i)
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE spectr1d

!*****************************************************************
      SUBROUTINE entrans(a,b,c,d,e,f,nmb,kin)
!-----------------------------------------------------------------
!
! Computes the energy (or cross-helicity) transfer 
! in Fourier space in 3D. The output is written to 
! a file by the first node.
!
! Output files contain:
! 'ktransfer.XXX.txt': k, Tv(k) (kinetic energy transfer function)
! 'mtransfer.XXX.txt': k, Tb(k) (magnetic energy transfer)
! 'jtransfer.XXX.txt': k, Tj(k) (Lorentz force work)
! 'kcrostran.XXX.txt': k, TC_v(k) (kinetic cross-helicity transfer)
! 'mcrostran.XXX.txt': k, TC_b(k) (magnetic cross-helicity transfer)
!
! Parameters
!     a  : field component in the x-direction
!     b  : field component in the y-direction
!     c  : field component in the z-direction
!     d  : nonlinear term in the x-direction
!     e  : nonlinear term in the y-direction
!     f  : nonlinear term in the z-direction
!     nmb: the extension used when writting the file
!     kin: =0 computes the magnetic energy transfer
!          =1 computes the kinetic energy transfer
!          =2 computes the Lorentz force work (energy transfer)
!          =3 computes the magnetic cross-helicity transfer
!          =4 computes the kinetic cross-helicity transfer
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nmax/2+1) :: Ek,Ektot
      DOUBLE PRECISION    :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: d,e,f
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,nmax/2+1
         Ek(i) = 0.0D0
      END DO
!
! Computes the kinetic energy transfer
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      IF (kin.ge.1) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(sqrt(kk2(k,j,1))/Dkk+.501)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = (real(a(k,j,1)*conjg(d(k,j,1)))+            &
                            real(b(k,j,1)*conjg(e(k,j,1)))+            &
                            real(c(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(sqrt(kk2(k,j,1))/Dkk+.501)
                     IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(d(k,j,i)))+       &
                                 real(b(k,j,i)*conjg(e(k,j,i)))+       &
                                 real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                 END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(sqrt(kk2(k,j,1))/Dkk+.501)
                     IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(d(k,j,i)))+       &
                                 real(b(k,j,i)*conjg(e(k,j,i)))+       &
                                 real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the magnetic energy transfer
!
      ELSE
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(sqrt(kk2(k,j,1))/Dkk+.501)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = kk2(k,j,1)*                                 &
                           (real(a(k,j,1)*conjg(d(k,j,1)))+            &
                            real(b(k,j,1)*conjg(e(k,j,1)))+            &
                            real(c(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(sqrt(kk2(k,j,1))/Dkk+.501)
                     IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                        tmq = 2*kk2(k,j,i)*                            &
                              (real(a(k,j,i)*conjg(d(k,j,i)))+         &
                               real(b(k,j,i)*conjg(e(k,j,i)))+         &
                               real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic 
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                 END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(sqrt(kk2(k,j,1))/Dkk+.501)
                     IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                        tmq = 2*kk2(k,j,i)*                            &
                              (real(a(k,j,i)*conjg(d(k,j,i)))+         &
                               real(b(k,j,i)*conjg(e(k,j,i)))+         &
                               real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Ek,Ektot,nmax/2+1,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF (kin.eq.0) THEN
            OPEN(1,file='mtransfer.' // nmb // '.txt')
         ELSEIF (kin.eq.1) THEN
            OPEN(1,file='ktransfer.' // nmb // '.txt')
         ELSEIF (kin.eq.2) THEN
            OPEN(1,file='jtransfer.' // nmb // '.txt')
         ELSEIF (kin.eq.3) THEN
            OPEN(1,file='mcrostran.' // nmb // '.txt')
         ELSEIF (kin.eq.4) THEN
            OPEN(1,file='kcrostran.' // nmb // '.txt')
         ENDIF
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)')  Dkk*i,Ektot(i)
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE entrans

!*****************************************************************
      SUBROUTINE heltrans(a,b,c,d,e,f,nmb,kin)
!-----------------------------------------------------------------
!
! Computes the helicity transfer in Fourier 
! space in 3D. The output is written to a 
! file by the first node.
!
! Output files contain:
! 'hktransfer.XXX.txt': k, TH_v(k) (transfer of kinetic helicity)
! 'hmtransfer.XXX.txt': k, TH_b(k) (transfer of magnetic helicity)
!
! Parameters
!     a  : field component in the x-direction (v or a)
!     b  : field component in the y-direction (v or a)
!     c  : field component in the z-direction (v or a)
!     d  : nonlinear term in the x-direction
!     e  : nonlinear term in the y-direction
!     f  : nonlinear term in the z-direction
!     nmb: the extension used when writting the file
!     kin: =0 computes the magnetic helicity transfer
!          =1 computes the kinetic helicity transfer
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nmax/2+1) :: Hk,Hktot
      DOUBLE PRECISION    :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)          :: c1,c2,c3
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Hk to zero
!
      DO i = 1,nmax/2+1
         Hk(i) = 0.0D0
      END DO
!
! Computes the helicity transfer
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      CALL rotor3(b,c,c1,1)
      CALL rotor3(a,c,c2,2)
      CALL rotor3(a,b,c3,3)
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,ny
            DO k = 1,nz
               kmn = int(sqrt(kk2(k,j,1))/Dkk+.501)
               IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                  tmq = (real(c1(k,j,1)*conjg(d(k,j,1)))+            &
                         real(c2(k,j,1)*conjg(e(k,j,1)))+            &
                         real(c3(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                  Hk(kmn) = Hk(kmn)+tmq
               ENDIF
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(sqrt(kk2(k,j,1))/Dkk+.501)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = 2*(real(c1(k,j,i)*conjg(d(k,j,i)))+       &
                              real(c2(k,j,i)*conjg(e(k,j,i)))+       &
                              real(c3(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                     Hk(kmn) = Hk(kmn)+tmq
                  ENDIF
              END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(sqrt(kk2(k,j,1))/Dkk+.501)
                  IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                     tmq = 2*(real(c1(k,j,i)*conjg(d(k,j,i)))+       &
                              real(c2(k,j,i)*conjg(e(k,j,i)))+       &
                              real(c3(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                     Hk(kmn) = Hk(kmn)+tmq
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Hk,Hktot,nmax/2+1,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF (kin.eq.0) THEN
            OPEN(1,file='hmtransfer.' // nmb // '.txt')
         ELSE
            OPEN(1,file='hktransfer.' // nmb // '.txt')
         ENDIF
         DO i=1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)')  Dkk*i,Hktot(i)
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE heltrans

!*****************************************************************
      SUBROUTINE pspectrum(a,fnout,svmax)
!-----------------------------------------------------------------
!
! Computes the 'power' spectrum of specified 
! scalar quantity. The output is written to a 
! file by the first node.
!
! Parameters
!     a    : input complex array
!     fnout: output file name 
!     svmax: number of spectral values to output
!
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nmax/2+1)                    :: Ek
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a
      INTEGER         , INTENT(IN)                           :: svmax
      INTEGER                                                :: i,imax
      CHARACTER(len=*), INTENT(IN)                           :: fnout
!
! Computes the energy and/or helicity spectra
      CALL pspectrumc(a,Ek)
!
! Exports the energy spectrum to a file
!
      IF (myrank.eq.0) THEN
         imax = min(svmax,nmax/2+1)
         OPEN(1,file=trim(fnout))
         DO i=1,imax
            WRITE(1,FMT='(E13.6,E23.15)')  Dkk*i,Ek(i)
         END DO
         CLOSE(1)
      ENDIF
!
      RETURN
      END SUBROUTINE pspectrum

!*****************************************************************
      SUBROUTINE pspectrumc(a,Ektot)
!-----------------------------------------------------------------
!
! Computes the 'power' spectrum of a, returning it
!
! Parameters
!     a    : input array
!     Ektot: output power spectrum
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION,              DIMENSION(nmax/2+1)      :: Ek
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(nmax/2+1)      :: Ektot
      DOUBLE PRECISION                                        :: tmq
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)              :: c1
      REAL(KIND=GP)       :: tmp
      INTEGER             :: i,j,k
      INTEGER             :: kmn

!
! Computes the power spectrum
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      DO i = 1,nmax/2+1
         Ek   (i) = 0.0D0
         Ektot(i) = 0.0D0
      END DO
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,ny
            DO k = 1,nz
               kmn = int(sqrt(kk2(k,j,1))/Dkk+.501)
                IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                   tmq = (abs(a(k,j,1))**2 )*tmp
!$omp atomic
                   Ek(kmn) = Ek(kmn)+tmq
                ENDIF
             END DO
          END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
          DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
             DO j = 1,ny
                DO k = 1,nz
                   kmn = int(sqrt(kk2(k,j,1))/Dkk+.501)
                   IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                      tmq = 2*(abs(a(k,j,i))**2 )*tmp
!$omp atomic
                      Ek(kmn) = Ek(kmn)+tmq
                   ENDIF
                END DO
             END DO
          END DO
        ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
          DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
             DO j = 1,ny
                DO k = 1,nz
                   kmn = int(sqrt(kk2(k,j,1))/Dkk+.501)
                   IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
                      tmq = 2*(abs(a(k,j,i))**2)*tmp
!$omp atomic
                      Ek(kmn) = Ek(kmn)+tmq
                   ENDIF
                END DO
             END DO
          END DO
       ENDIF
!
! Computes the reduction between nodes

       CALL MPI_ALLREDUCE(Ek,Ektot,nmax/2+1,MPI_DOUBLE_PRECISION, &
                      MPI_SUM,MPI_COMM_WORLD,ierr)
!
       RETURN
       END SUBROUTINE pspectrumc

!*****************************************************************
      SUBROUTINE phaseshift(a,cdump)
!-----------------------------------------------------------------
!
! Shifts an array 'a' in complex space by a phase. This subroutine
! computes the shifted array while taking care of the planes in
! Fourier space with kx=0, ky=0, and ky=0, which must be dealt
! with care to keep the Fourier transform of 'a' real.
!
! Parameters
!     a    : input/output array
!     cdump: e(i.phase), where 'phase' is the phase shift
!
      USE fprecision
      USE mpivars
      USE grid
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(IN) :: cdump
      COMPLEX(KIND=GP)             :: jdump      
      INTEGER                      :: i,j,k

      jdump = conjg(cdump)
      IF (ista.eq.1) THEN
!$omp parallel do
         DO j = 2,ny/2+1
            a(1,j,1) = a(1,j,1)*cdump
            a(1,ny-j+2,1) = a(1,ny-j+2,1)*jdump
         END DO
!$omp parallel do
         DO k = 2,nz/2+1
            a(k,1,1) = a(k,1,1)*cdump
            a(nz-k+2,1,1) = a(nz-k+2,1,1)*jdump
         END DO
!$omp parallel do private (k)
         DO j = 2,ny
            DO k = 2,nz/2+1
               a(k,j,1) = a(k,j,1)*cdump
               a(nz-k+2,ny-j+2,1) = a(nz-k+2,ny-j+2,1)*jdump
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  a(k,j,i) = a(k,j,i)*cdump
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  a(k,j,i) = a(k,j,i)*cdump
               END DO
            END DO
         END DO
      ENDIF

      RETURN
      END SUBROUTINE phaseshift
