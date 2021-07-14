!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Extra subroutines to compute spatial derivatives and 
! nonlinear terms in compressible HD and MHD equations in 3D using 
! a pseudo-spectral method. You should use the FFTPLANS and 
! MPIVARS modules (see the file 'fftp_mod.f90') in each 
! program that calls any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2015 Pablo Dmitruk.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: pdmitruk@df.uba.ar 
!=================================================================

!*****************************************************************
      SUBROUTINE divide(d,a,b,c)
!-----------------------------------------------------------------
!
! Computes the division of vector field A=(a,b,c) 
! by scalar field 'd' in real space.  
!
! Parameters
!     a: input/output matrix with A_x/(A_x/d) (in Fourier space)
!     b: input/output matrix with A_y/(A_y/d) (in Fourier space)
!     c: input/output matrix with A_z/(A_z/d) (in Fourier space)
!     d: input matrix with scalar field
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: c
      COMPLEX(KIND=GP), INTENT(IN),    DIMENSION(nz,ny,ista:iend) :: d
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: x
      REAL(KIND=GP),    DIMENSION(nx,ny,ksta:kend) :: r1,r2
      REAL(KIND=GP),    DIMENSION(nx,ny,ksta:kend) :: r3,r4
      INTEGER :: i,j,k

      CALL fftp3d_complex_to_real(plancr,a,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,b,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c,r3,MPI_COMM_WORLD)

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               x(k,j,i) = d(k,j,i)
            END DO
         END DO
      END DO

      CALL fftp3d_complex_to_real(plancr,x,r4,MPI_COMM_WORLD)

!$omp parallel do if (iend-ista.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (iend-ista.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r1(i,j,k) = r1(i,j,k)/r4(i,j,k)
               r2(i,j,k) = r2(i,j,k)/r4(i,j,k)
               r3(i,j,k) = r3(i,j,k)/r4(i,j,k)
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,r1,a,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r2,b,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r3,c,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE divide

!*****************************************************************
      SUBROUTINE gradpress(cp1,gam1,d,a,b,c,e,f,g)
!-----------------------------------------------------------------
!
! Computes the gradient of the pressure = 0.5(vel*vel + cp1*rho^gam1) 
! with gam1 = gamma - 1
!
! Parameters
!     a  : input matrix with v_x (in Fourier space)
!     b  : input matrix with v_y (in Fourier space)
!     c  : input matrix with v_z (in Fourier space)
!     d  : input matrix with density (in Fourier space)
!     e  : output matrix with grad(press)_x (in Fourier space)
!     f  : output matrix with grad(press)_y (in Fourier space)
!     g  : output matrix with grad(press)_z (in Fourier space)
!
      USE fprecision
      USE kes
      USE grid
      USE commtypes
      USE mpivars
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(nz,ny,ista:iend) :: c,d
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: e,f
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: g
      REAL(KIND=GP),    INTENT(IN)                 :: cp1, gam1
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: h
      REAL(KIND=GP),    DIMENSION(nx,ny,ksta:kend) :: r1,r2
      REAL(KIND=GP),    DIMENSION(nx,ny,ksta:kend) :: r3,r4
      REAL(KIND=GP) :: tmp
      INTEGER       :: i,j,k

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               h(k,j,i) = d(k,j,i)
               e(k,j,i) = a(k,j,i)
               f(k,j,i) = b(k,j,i)
               g(k,j,i) = c(k,j,i)
            END DO
         END DO
      END DO

      CALL fftp3d_complex_to_real(plancr,e,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,f,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,g,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,h,r4,MPI_COMM_WORLD)

      tmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (iend-ista.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = 0.5_GP*( r1(i,j,k)*r1(i,j,k)*tmp*tmp + &
                                    r2(i,j,k)*r2(i,j,k)*tmp*tmp + &
                                    r3(i,j,k)*r3(i,j,k)*tmp*tmp + &
                           cp1*((r4(i,j,k)*tmp)**gam1) )
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,r4,h,MPI_COMM_WORLD)
      CALL derivk3(h,e,1)
      CALL derivk3(h,f,2)
      CALL derivk3(h,g,3)

      RETURN
      END SUBROUTINE gradpress

!*****************************************************************
      SUBROUTINE gradpstate(cp1,gam1,d,e,f,g)
!-----------------------------------------------------------------
!
! Computes the gradient of the pressure resulting only from
! the equation of state of the gas = 0.5*cp1*rho^gam1, with
! gam1 = gamma - 1
!
! Parameters
!     d  : input matrix with density (in Fourier space)
!     e  : output matrix with grad(press)_x (in Fourier space)
!     f  : output matrix with grad(press)_y (in Fourier space)
!     g  : output matrix with grad(press)_z (in Fourier space)
!
      USE fprecision
      USE kes
      USE grid
      USE commtypes
      USE mpivars
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(nz,ny,ista:iend) :: d
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: e,f,g
      REAL(KIND=GP),    INTENT(IN)                 :: cp1, gam1
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: h
      REAL(KIND=GP),    DIMENSION(nx,ny,ksta:kend) :: r4
      REAL(KIND=GP) :: tmp
      INTEGER       :: i,j,k

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               h(k,j,i) = d(k,j,i)
            END DO
         END DO
      END DO

      CALL fftp3d_complex_to_real(plancr,h,r4,MPI_COMM_WORLD)

      tmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (iend-ista.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = .5_GP*cp1*(r4(i,j,k)*tmp)**gam1
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,r4,h,MPI_COMM_WORLD)
      CALL derivk3(h,e,1)
      CALL derivk3(h,f,2)
      CALL derivk3(h,g,3)

      RETURN
      END SUBROUTINE gradpstate

!*****************************************************************
      SUBROUTINE divrhov(d,a,b,c,e)
!-----------------------------------------------------------------
!
! Computes the divergence of the product of scalar 'd' by
! vector A = (a,b,c) 
!
! Parameters
!     a  : input matrix with v_x (in Fourier space)
!     b  : input matrix with v_y (in Fourier space)
!     c  : input matrix with v_z (in Fourier space) [A = (a,b,c)]
!     d  : input matrix with density (in Fourier space)
!     e  : output matrix with div(d.A) (in Fourier space)
!
      USE fprecision
      USE kes
      USE grid
      USE commtypes
      USE mpivars
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(nz,ny,ista:iend) :: c,d
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: e
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: f,g,h
      REAL(KIND=GP),    DIMENSION(nx,ny,ksta:kend) :: r1,r2
      REAL(KIND=GP),    DIMENSION(nx,ny,ksta:kend) :: r3,r4
      REAL(KIND=GP) :: tmp
      INTEGER       :: i,j,k

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               e(k,j,i) = d(k,j,i)
               f(k,j,i) = a(k,j,i)
               g(k,j,i) = b(k,j,i)
               h(k,j,i) = c(k,j,i)
            END DO
         END DO
      END DO

      CALL fftp3d_complex_to_real(plancr,f,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,g,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,h,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,e,r4,MPI_COMM_WORLD)

      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!$omp parallel do if (iend-ista.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (iend-ista.lt.nth) private (i)
         DO j = 1,ny
            DO i = 1,nx
               r1(i,j,k) = r4(i,j,k)*r1(i,j,k)*tmp
               r2(i,j,k) = r4(i,j,k)*r2(i,j,k)*tmp
               r3(i,j,k) = r4(i,j,k)*r3(i,j,k)*tmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,r1,f,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r2,g,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r3,h,MPI_COMM_WORLD)
      CALL derivk3(f,e,1)
      CALL derivk3(g,f,2)
      CALL derivk3(h,g,3)

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               e(k,j,i) = e(k,j,i) + f(k,j,i) + g(k,j,i)
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE divrhov

!*****************************************************************
      SUBROUTINE vdiss(nu,nu2,a,b,c)
!-----------------------------------------------------------------
!
! Computes the kinetic dissipation term 
! nu*del^2(vel) + nu2*grad(div(vel))
!
! Parameters
!     a  : input/output matrix with v_x/(diss(v)_x) (in Fourier space)
!     b  : input/output matrix with v_y/(diss(v)_y) (in Fourier space)
!     c  : input/output matrix with v_z/(diss(v)_z) (in Fourier space)
!
      USE fprecision
      USE kes
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: b
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: c
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: d,e,f,g
      REAL(KIND=GP),    INTENT(IN)                 :: nu,nu2
      INTEGER                                      :: i,j,k

                                          ! div(vel)
      CALL derivk3(a,e,1)
      CALL derivk3(b,f,2)
      CALL derivk3(c,g,3)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               d(k,j,i) = e(k,j,i)+f(k,j,i)+g(k,j,i)
            END DO
         END DO
      END DO
      CALL derivk3(d,e,1)
      CALL derivk3(d,f,2)
      CALL derivk3(d,g,3)
                                          ! del^2(vel)
      CALL laplak3(a,a)
      CALL laplak3(b,b)
      CALL laplak3(c,c)
      
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               a(k,j,i) = nu*a(k,j,i)+nu2*e(k,j,i)
               b(k,j,i) = nu*b(k,j,i)+nu2*f(k,j,i)
               c(k,j,i) = nu*c(k,j,i)+nu2*g(k,j,i)
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE vdiss

!*****************************************************************
      SUBROUTINE energycompr(gam1,cp1, d,a,b,c,t,dt)
!-----------------------------------------------------------------
!
! Computes the kinetic and internal energy for compressible runs, 
! including the mass density 
!
! Output file contains:
! 'compr_ener.txt': time, kinetic energy, internal energy
!
! Parameters
!     a  : input matrix with v_x (in Fourier space)
!     b  : input matrix with v_y (in Fourier space)
!     c  : input matrix with v_z (in Fourier space) [A = (a,b,c)]
!     d  : input matrix with density (in Fourier space)
!     t  : number of time steps made
!     dt : time step
!     gam1 : gamma-1 constant (adiabatic constant)
!
      USE fprecision
      USE kes
      USE grid
      USE commtypes
      USE mpivars
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: b
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: d
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: e,f,g,h
      REAL(KIND=GP),    DIMENSION(nx,ny,ksta:kend) :: r1,r2
      REAL(KIND=GP),    DIMENSION(nx,ny,ksta:kend) :: r3,r4
      REAL(KIND=GP),    INTENT(IN)  :: gam1,cp1,dt
      INTEGER,          INTENT(IN)  :: t
      REAL(KIND=GP)                 :: tmp, tmp1, gam0
      DOUBLE PRECISION              :: tot_ekin
      DOUBLE PRECISION              :: loc_ekin
      DOUBLE PRECISION              :: tot_eint
      DOUBLE PRECISION              :: loc_eint
      INTEGER                       :: i,j,k

      tot_ekin = 0.0D0
      tot_eint = 0.0D0

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               e(k,j,i) = d(k,j,i)
               f(k,j,i) = a(k,j,i)
               g(k,j,i) = b(k,j,i)
               h(k,j,i) = c(k,j,i)
            END DO
         END DO
      END DO

      CALL fftp3d_complex_to_real(plancr,f,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,g,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,h,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,e,r4,MPI_COMM_WORLD)

      loc_ekin = 0.0D0
      loc_eint = 0.0D0
      tmp  = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**3
      tmp1 = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
      gam0 = gam1 + 1.0_GP
!$omp parallel do if (iend-ista.ge.nth) private (j,i) reduction(+:loc_ekin,loc_eint)
      DO k = ksta,kend
!$omp parallel do if (iend-ista.lt.nth) private (i) reduction(+:loc_ekin,loc_eint)
         DO j = 1,ny
            DO i = 1,nx
               loc_ekin = loc_ekin + r4(i,j,k) * (r1(i,j,k)*r1(i,j,k)   + &
                                                  r2(i,j,k)*r2(i,j,k)   + &
                                                  r3(i,j,k)*r3(i,j,k) ) * tmp
               loc_eint = loc_eint + (r4(i,j,k)*tmp1)**gam0
            END DO
         END DO
      END DO

      loc_ekin = loc_ekin*tmp1
      loc_eint = loc_eint*tmp1

      CALL MPI_REDUCE(loc_ekin,tot_ekin,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(loc_eint,tot_eint,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      tot_eint = tot_eint*cp1/gam0 

      IF (myrank.eq.0) THEN
         OPEN(1,file='compr_energy.txt',position='append')
         WRITE(1,*) (t-1)*dt,tot_ekin,tot_eint
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE energycompr
