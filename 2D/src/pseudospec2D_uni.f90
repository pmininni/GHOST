!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines for computing spatial derivatives and nonlinear 
! terms in incompressible HD, MHD and Hall-MHD equations with 
! a uniform magnetic field in 2D using a pseudo-spectral method. 
! You should use the FFTPLANS and MPIVARS modules (see the file 
! 'fftp2D_mod.f90') in each program that call any of the 
! subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!
! 2003 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!=================================================================

!*****************************************************************
      SUBROUTINE poissonb0(a,b,c,b0)
!-----------------------------------------------------------------
!
! Poisson bracket of the scalar fields A and B 
! in real space.
!
! Parameters
!     a : input matrix
!     b : input matrix
!     c : Poisson bracket {A,B} [output]
!     b0: Amplitude of the uniform field in y
!
      USE mpivars
      USE grid
      USE fft
      IMPLICIT NONE

      COMPLEX, DIMENSION(n,ista:iend) :: a,b,c
      COMPLEX, DIMENSION(n,ista:iend) :: c1,c2
      REAL, DIMENSION(n,jsta:jend)    :: r1,r2,r3
      REAL    :: b0
      INTEGER :: i,j

!
! Computes dA/dx.dB/dy
!
      CALL derivk2(a,c1,1)
      CALL derivk2(b,c2,2)
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,n
            r3(i,j) = r1(i,j)*r2(i,j)
         END DO
      END DO
!
! Computes (dA/dy+b0).dB/dx
!
      CALL derivk2(a,c1,2)
      CALL derivk2(b,c2,1)
      c2(1,1) = -b0*float(n)**2
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,n
            r3(i,j) = (r3(i,j)-r1(i,j)*r2(i,j))/float(n)**4
         END DO
      END DO

      CALL fftp2d_real_to_complex(planrc,r3,c,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE poissonb0
