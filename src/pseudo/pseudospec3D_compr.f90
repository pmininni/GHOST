!=================================================================
! PSEUDOSPECTRAL modules
!
! CONTAINS:
!      MODULE pseudospec_compressible
!      MODULE pseudospec_compr
!
! Extra subroutines to compute spatial derivatives and
! nonlinear terms in compressible HD and MHD equations in 3D using
! a pseudo-spectral method. You should use the FFTPLANS and
! MPIVARS modules (see the file 'fftp_mod.f90') in each
! program that calls any of the subroutines in this file.
!
! The module pseudospec_compressible has the routines used in the
! evolution equations (kernels with the dual host/device
! directives, temporaries from the workspace pool), and the
! module pseudospec_compr the diagnostics (computed in the host,
! with host-only temporaries). See pseudospec3D_hd.f90 for the
! methodology.
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

MODULE pseudospec_compressible
   USE pseudospec_fluid
   USE class_GWorkspace3D, ONLY: gws
   USE gdevice, ONLY: gdev_active
   CONTAINS

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
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN),    DIMENSION(nz,ny,ista:iend) :: d
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: x
      REAL(KIND=GP)   , POINTER, DIMENSION(:,:,:) :: r1,r2,r3,r4
      INTEGER :: i,j,k
      LOGICAL :: bret

      CALL gws%get_complex_tmp(x,bret)
      CALL gws%get_real_tmp(r1,bret)
      CALL gws%get_real_tmp(r2,bret)
      CALL gws%get_real_tmp(r3,bret)
      CALL gws%get_real_tmp(r4,bret)

      CALL fftp3d_complex_to_real(plancr,a,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,b,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,c,r3,MPI_COMM_WORLD)
      CALL copy3(d,x)
      CALL fftp3d_complex_to_real(plancr,x,r4,MPI_COMM_WORLD)

#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
      DO k = ksta,kend
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

      CALL gws%free_real_tmp(r4)
      CALL gws%free_real_tmp(r3)
      CALL gws%free_real_tmp(r2)
      CALL gws%free_real_tmp(r1)
      CALL gws%free_complex_tmp(x)

      RETURN
      END SUBROUTINE divide

!*****************************************************************
      SUBROUTINE gradpressi(gam1,e,dpx,dpy,dpz)
!-----------------------------------------------------------------
!
! Computes the gradient of the thermo. pressure assuming
! ideal equation of state, based on a polytropic law, s.t.
!    p = (gamma - 1) * e
! where e is internal energy density.
!
! Parameters
!     gam1: gamma - 1
!     e   : input matrix with int. energy density (in Fourier space)
!     dpx : output matrix with grad(press)_x (in Fourier space)
!     dpy : output matrix with grad(press)_y (in Fourier space)
!     dpz : output matrix with grad(press)_z (in Fourier space)
!
      USE fprecision
      USE kes
      USE grid
      USE commtypes
      USE mpivars
      USE fft
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(nz,ny,ista:iend) :: e
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: dpx,dpy,dpz
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: t
      REAL(KIND=GP)   , INTENT(IN)                :: gam1
      INTEGER       :: i,j,k
      LOGICAL       :: bret

      CALL gws%get_complex_tmp(t,bret)
      CALL copy3(e,t)
      CALL scal3(t,gam1)
      CALL derivk3(t,dpx,1)
      CALL derivk3(t,dpy,2)
      CALL derivk3(t,dpz,3)
      CALL gws%free_complex_tmp(t)

      RETURN
      END SUBROUTINE gradpressi

!*****************************************************************
      SUBROUTINE gradpress(cp1,gam1,d,a,b,c,e,f,g)
!-----------------------------------------------------------------
!
! Computes the gradient of the pressure = 0.5(vel*vel + cp1*rho^gam1)
! with gam1 = gamma - 1 (the sum of the kinetic energy per unit
! mass and of the enthalpy of a polytropic gas)
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
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(nz,ny,ista:iend) :: c,d
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: e,f
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: g
      REAL(KIND=GP)   , INTENT(IN)                :: cp1, gam1
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: h
      REAL(KIND=GP)   , POINTER, DIMENSION(:,:,:) :: r1,r2,r3,r4
      REAL(KIND=GP) :: tmp
      INTEGER       :: i,j,k
      LOGICAL       :: bret

      CALL gws%get_complex_tmp(h,bret)
      CALL gws%get_real_tmp(r1,bret)
      CALL gws%get_real_tmp(r2,bret)
      CALL gws%get_real_tmp(r3,bret)
      CALL gws%get_real_tmp(r4,bret)

      CALL copy3(a,e)
      CALL copy3(b,f)
      CALL copy3(c,g)
      CALL copy3(d,h)
      CALL fftp3d_complex_to_real(plancr,e,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,f,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,g,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,h,r4,MPI_COMM_WORLD)

      tmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
      DO k = ksta,kend
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

      CALL gws%free_real_tmp(r4)
      CALL gws%free_real_tmp(r3)
      CALL gws%free_real_tmp(r2)
      CALL gws%free_real_tmp(r1)
      CALL gws%free_complex_tmp(h)

      RETURN
      END SUBROUTINE gradpress

!*****************************************************************
      SUBROUTINE gradpstate(cp1,gam1,d,e,f,g)
!-----------------------------------------------------------------
!
! Computes the gradient of the pressure resulting only from
! the equation of state of the gas = 0.5*cp1*rho^gam1, with
! gam1 = gamma - 1 (the enthalpy of a polytropic gas)
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
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(nz,ny,ista:iend) :: d
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: e,f,g
      REAL(KIND=GP),    INTENT(IN)                :: cp1, gam1
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: h
      REAL(KIND=GP)   , POINTER, DIMENSION(:,:,:) :: r4
      REAL(KIND=GP) :: tmp
      INTEGER       :: i,j,k
      LOGICAL       :: bret

      CALL gws%get_complex_tmp(h,bret)
      CALL gws%get_real_tmp(r4,bret)

      CALL copy3(d,h)
      CALL fftp3d_complex_to_real(plancr,h,r4,MPI_COMM_WORLD)

      tmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
      DO k = ksta,kend
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

      CALL gws%free_real_tmp(r4)
      CALL gws%free_complex_tmp(h)

      RETURN
      END SUBROUTINE gradpstate

!*****************************************************************
      SUBROUTINE divrhov(d,a,b,c,dodealias,e)
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
!     dodealias:  flag (0, 1) to do dealiasing
!     e  : output matrix with div(d.A) (in Fourier space)
!
      USE fprecision
      USE kes
      USE grid
      USE commtypes
      USE mpivars
      USE fft
      USE ali
      IMPLICIT NONE

      INTEGER         , INTENT (IN)                             :: dodealias
      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(nz,ny,ista:iend) :: c,d
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: e
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: f,g,h
      REAL(KIND=GP)   , POINTER, DIMENSION(:,:,:) :: r1,r2,r3,r4
      REAL(KIND=GP) :: tmp
      INTEGER       :: i,j,k
      LOGICAL       :: bret

      CALL gws%get_complex_tmp(f,bret)
      CALL gws%get_complex_tmp(g,bret)
      CALL gws%get_complex_tmp(h,bret)
      CALL gws%get_real_tmp(r1,bret)
      CALL gws%get_real_tmp(r2,bret)
      CALL gws%get_real_tmp(r3,bret)
      CALL gws%get_real_tmp(r4,bret)

      CALL copy3(d,e)
      CALL copy3(a,f)
      CALL copy3(b,g)
      CALL copy3(c,h)

      ! Dealiases the fields:
      IF ( dodealias .gt. 0 ) THEN
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
                  IF (kn2(k,j,i).gt.kmax) THEN
                     e(k,j,i) = 0.0_GP
                     f(k,j,i) = 0.0_GP
                     g(k,j,i) = 0.0_GP
                     h(k,j,i) = 0.0_GP
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF

      CALL fftp3d_complex_to_real(plancr,f,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,g,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,h,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,e,r4,MPI_COMM_WORLD)

      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
      DO k = ksta,kend
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

#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               e(k,j,i) = e(k,j,i) + f(k,j,i) + g(k,j,i)
            END DO
         END DO
      END DO

      ! Dealiases the result:
      IF ( dodealias .gt. 0 ) THEN
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
                  IF (kn2(k,j,i).gt.kmax) THEN
                     e(k,j,i) = 0.0_GP
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF

      CALL gws%free_real_tmp(r4)
      CALL gws%free_real_tmp(r3)
      CALL gws%free_real_tmp(r2)
      CALL gws%free_real_tmp(r1)
      CALL gws%free_complex_tmp(h)
      CALL gws%free_complex_tmp(g)
      CALL gws%free_complex_tmp(f)

      RETURN
      END SUBROUTINE divrhov

!*****************************************************************
      SUBROUTINE vdiss(nu,nu2,a,b,c,d,e,f)
!-----------------------------------------------------------------
!
! Computes the kinetic dissipation term
! nu*del^2(vel) + nu2*grad(div(vel))
!
! Parameters
!     nu : kinematic viscosity
!     nu2: second (bulk) viscosity for the divergence term
!     a  : input matrix with v_x (in Fourier space)
!     b  : input matrix with v_y (in Fourier space)
!     c  : input matrix with v_z (in Fourier space)
!     d  : output matrix with diss(v)_x (in Fourier space)
!     e  : output matrix with diss(v)_y (in Fourier space)
!     f  : output matrix with diss(v)_z (in Fourier space)
!
      USE fprecision
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: g,h
      REAL(KIND=GP)   , INTENT(IN)                :: nu,nu2
      INTEGER                                     :: i,j,k
      LOGICAL                                     :: bret

      CALL gws%get_complex_tmp(g,bret)
      CALL gws%get_complex_tmp(h,bret)
                                          ! div(vel)
      CALL derivk3(a,d,1)
      CALL derivk3(b,e,2)
      CALL derivk3(c,f,3)
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               g(k,j,i) = d(k,j,i)+e(k,j,i)+f(k,j,i)
            END DO
         END DO
      END DO
                                         ! nu del^2(vel) + nu2 grad(div(vel))
      CALL derivk3(g,h,1)
      CALL laplak3(a,d)
      CALL saxpby_c(d,d,nu,h,nu2)
      CALL derivk3(g,h,2)
      CALL laplak3(b,e)
      CALL saxpby_c(e,e,nu,h,nu2)
      CALL derivk3(g,h,3)
      CALL laplak3(c,f)
      CALL saxpby_c(f,f,nu,h,nu2)

      CALL gws%free_complex_tmp(h)
      CALL gws%free_complex_tmp(g)

      RETURN
      END SUBROUTINE vdiss

!*****************************************************************
      SUBROUTINE pdVwork(gam1,e,a,b,c,dodealias,pdV)
!-----------------------------------------------------------------
!
! Computes p.Div v term
!
! Parameters
!     gam1: gamma - 1
!     e   : input matrix with int. energy density (in Fourier space)
!     a   : input matrix with v_x (in Fourier space)
!     b   : input matrix with v_y (in Fourier space)
!     c   : input matrix with v_z (in Fourier space) [A = (a,b,c)]
!     dodealias: flag (0, 1) to do dealiasing
!     pdV : result
!
      USE fprecision
      USE kes
      USE grid
      USE commtypes
      USE mpivars
      USE fft
      USE ali
      IMPLICIT NONE

      INTEGER         , INTENT (IN)                             :: dodealias
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: b
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: c
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: e
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: pdV
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: t1,t2,t3,t4
      REAL(KIND=GP)   , INTENT (IN)               :: gam1
      REAL(KIND=GP)   , POINTER, DIMENSION(:,:,:) :: r1,r2,r3,r4
      REAL(KIND=GP)                 :: tmp
      INTEGER                       :: i,j,k
      LOGICAL                       :: bret

      CALL gws%get_complex_tmp(t1,bret)
      CALL gws%get_complex_tmp(t2,bret)
      CALL gws%get_complex_tmp(t3,bret)
      CALL gws%get_complex_tmp(t4,bret)
      CALL gws%get_real_tmp(r1,bret)
      CALL gws%get_real_tmp(r2,bret)
      CALL gws%get_real_tmp(r3,bret)
      CALL gws%get_real_tmp(r4,bret)

      ! Take divergence terms:
      CALL derivk3(a,t1,1)
      CALL derivk3(b,t2,2)
      CALL derivk3(c,t3,3)
      CALL copy3(e,t4)

      ! Dealiases the fields:
      IF ( dodealias .gt. 0 ) THEN
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
                  IF (kn2(k,j,i).gt.kmax) THEN
                     t1(k,j,i) = 0.0_GP
                     t2(k,j,i) = 0.0_GP
                     t3(k,j,i) = 0.0_GP
                     t4(k,j,i) = 0.0_GP
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF

      CALL fftp3d_complex_to_real(plancr,t1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,t2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,t3,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,t4,r4,MPI_COMM_WORLD)

      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               r1(i,j,k) = gam1*r4(i,j,k)*(r1(i,j,k)+r2(i,j,k)+r3(i,j,k))*tmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,r1,pdV,MPI_COMM_WORLD)

      ! Dealiases the result:
      IF ( dodealias .gt. 0 ) THEN
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
                  IF (kn2(k,j,i).gt.kmax) THEN
                     pdV(k,j,i) = 0.0_GP
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF

      CALL gws%free_real_tmp(r4)
      CALL gws%free_real_tmp(r3)
      CALL gws%free_real_tmp(r2)
      CALL gws%free_real_tmp(r1)
      CALL gws%free_complex_tmp(t4)
      CALL gws%free_complex_tmp(t3)
      CALL gws%free_complex_tmp(t2)
      CALL gws%free_complex_tmp(t1)

      RETURN
      END SUBROUTINE pdVwork

!*****************************************************************
      SUBROUTINE mom2vel(rho,sx,sy,sz,dodealias,vx,vy,vz)
!-----------------------------------------------------------------
!
! Computes velocity from momentum
!
! Parameters
!     rho      : density
!     sx,sy, sz: momentum components
!     dodealias: flag (0, 1) to do dealiasing
!     vx,vy,vz : velocity components
!
      USE fprecision
      USE kes
      USE grid
      USE commtypes
      USE mpivars
      USE fft
      USE ali
      IMPLICIT NONE

      INTEGER         , INTENT (IN)                             :: dodealias
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: rho
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: sx,sy,sz
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: vx,vy,vz
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: t4
      REAL(KIND=GP)   , POINTER, DIMENSION(:,:,:) :: r1,r2,r3,r4
      INTEGER                       :: i,j,k
      LOGICAL                       :: bret

      CALL gws%get_complex_tmp(t4,bret)
      CALL gws%get_real_tmp(r1,bret)
      CALL gws%get_real_tmp(r2,bret)
      CALL gws%get_real_tmp(r3,bret)
      CALL gws%get_real_tmp(r4,bret)

      CALL copy3(sx ,vx)
      CALL copy3(sy ,vy)
      CALL copy3(sz ,vz)
      CALL copy3(rho,t4)

      ! Dealiases the fields:
      IF ( dodealias .gt. 0 ) THEN
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
                  IF (kn2(k,j,i).gt.kmax) THEN
                     vx(k,j,i) = 0.0_GP
                     vy(k,j,i) = 0.0_GP
                     vz(k,j,i) = 0.0_GP
                     t4(k,j,i) = 0.0_GP
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF

      CALL fftp3d_complex_to_real(plancr,vx,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,vy,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,vz,r3,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,t4,r4,MPI_COMM_WORLD)

#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               r1(i,j,k) = r1(i,j,k)/r4(i,j,k)
               r2(i,j,k) = r2(i,j,k)/r4(i,j,k)
               r3(i,j,k) = r3(i,j,k)/r4(i,j,k)
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,r1,vx,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r2,vy,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r3,vz,MPI_COMM_WORLD)

      ! Dealiases the result:
      IF ( dodealias .gt. 0 ) THEN
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
                  IF (kn2(k,j,i).gt.kmax) THEN
                     vx(k,j,i) = 0.0_GP
                     vy(k,j,i) = 0.0_GP
                     vz(k,j,i) = 0.0_GP
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF

      CALL gws%free_real_tmp(r4)
      CALL gws%free_real_tmp(r3)
      CALL gws%free_real_tmp(r2)
      CALL gws%free_real_tmp(r1)
      CALL gws%free_complex_tmp(t4)

      RETURN
      END SUBROUTINE mom2vel

!*****************************************************************
      SUBROUTINE viscHeatRayleigh(a,b,c,phi)
!-----------------------------------------------------------------
!
! Computes viscous heat term (Rayleigh form):
!     phi = tau_ij dv^i/dx^j
! where
!     tau_ij = 2 mu S_ij - 2/3 mu Div v delta_ij
! and
!     S_ij   = 1/2 ( v^i,j + v^j,i ) is strain rate.
! Then
!     phi = mu (v^j,i + v^i,j) - 2/3 mu (Div v)^2.
!
! Actually, the kernel phi/mu is returned, and user
! should multiply this term by mu
!
! Parameters
!     a   : input matrix with v_x (in Fourier space)
!     b   : input matrix with v_y (in Fourier space)
!     c   : input matrix with v_z (in Fourier space)
!     phi : result
!
      USE fprecision
      USE kes
      USE grid
      USE commtypes
      USE mpivars
      USE fft
      USE pseudospec_strain
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: b
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: phi
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: t1,t2,t3,t4,t5,t6,t7
      REAL(KIND=GP)   , POINTER, DIMENSION(:,:,:) :: r1,r2,r3,r4
      REAL(KIND=GP)                 :: tmp
      INTEGER                       :: i,j,k
      INTEGER                       :: btrunc,bnorm
      LOGICAL                       :: bret

      CALL gws%get_complex_tmp(t1,bret)
      CALL gws%get_complex_tmp(t2,bret)
      CALL gws%get_complex_tmp(t3,bret)
      CALL gws%get_complex_tmp(t4,bret)
      CALL gws%get_complex_tmp(t5,bret)
      CALL gws%get_complex_tmp(t6,bret)
      CALL gws%get_complex_tmp(t7,bret)
      CALL gws%get_real_tmp(r1,bret)
      CALL gws%get_real_tmp(r2,bret)
      CALL gws%get_real_tmp(r3,bret)
      CALL gws%get_real_tmp(r4,bret)

      ! Strain takes INOUT velocities: work on copies
      CALL copy3(a,t5)
      CALL copy3(b,t6)
      CALL copy3(c,t7)

      btrunc = 0
      bnorm  = 0
      ! Find diagonal strain rate components with normalization:
      CALL Strain(t5,t6,t7,1,1,btrunc,0.0_GP,0.0_GP,bnorm,t4,t1) ! S11
      CALL Strain(t5,t6,t7,2,2,btrunc,0.0_GP,0.0_GP,bnorm,t4,t2) ! S22
      CALL Strain(t5,t6,t7,3,3,btrunc,0.0_GP,0.0_GP,bnorm,t4,t3) ! S33

      CALL fftp3d_complex_to_real(plancr,t1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,t2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,t3,r3,MPI_COMM_WORLD)

      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = 2.0_GP*( r1(i,j,k)*r1(i,j,k) &
                                  + r2(i,j,k)*r2(i,j,k) &
                                  + r3(i,j,k)*r3(i,j,k) )*tmp
            END DO
         END DO
      END DO

      ! Find off-diagonal strain rate components with normalization:
      CALL Strain(t5,t6,t7,1,2,btrunc,0.0_GP,0.0_GP,bnorm,t4,t1) ! S12
      CALL Strain(t5,t6,t7,1,3,btrunc,0.0_GP,0.0_GP,bnorm,t4,t2) ! S13
      CALL Strain(t5,t6,t7,2,3,btrunc,0.0_GP,0.0_GP,bnorm,t4,t3) ! S23

      CALL fftp3d_complex_to_real(plancr,t1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,t2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,t3,r3,MPI_COMM_WORLD)

#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = r4(i,j,k) + 4.0_GP*( r1(i,j,k)*r1(i,j,k) &
                                              + r2(i,j,k)*r2(i,j,k) &
                                              + r3(i,j,k)*r3(i,j,k) )*tmp
            END DO
         END DO
      END DO

      ! Subtract dilatation term. First, compute it:
      CALL derivk3(a,t1,1)
      CALL derivk3(b,t2,2)
      CALL derivk3(c,t3,3)
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               t3(k,j,i) = t3(k,j,i) + t2(k,j,i) + t1(k,j,i)
            END DO
         END DO
      END DO

      CALL fftp3d_complex_to_real(plancr,t3,r3,MPI_COMM_WORLD)

#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               r4(i,j,k) = r4(i,j,k) - (2.0_GP/3.0_GP)*r3(i,j,k)*r3(i,j,k)*tmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,r4,phi,MPI_COMM_WORLD)

      CALL gws%free_real_tmp(r4)
      CALL gws%free_real_tmp(r3)
      CALL gws%free_real_tmp(r2)
      CALL gws%free_real_tmp(r1)
      CALL gws%free_complex_tmp(t7)
      CALL gws%free_complex_tmp(t6)
      CALL gws%free_complex_tmp(t5)
      CALL gws%free_complex_tmp(t4)
      CALL gws%free_complex_tmp(t3)
      CALL gws%free_complex_tmp(t2)
      CALL gws%free_complex_tmp(t1)

      RETURN
      END SUBROUTINE viscHeatRayleigh

END MODULE pseudospec_compressible


!=================================================================
! Diagnostics for compressible flows (computed in the host)
!=================================================================
MODULE pseudospec_compr
   USE pseudospec_fluid
   USE class_GWorkspace3D, ONLY: gws
   CONTAINS

!*****************************************************************
      SUBROUTINE energycompr(gam1,cp1,d,a,b,c,t,dt,path)
!-----------------------------------------------------------------
!
! Computes the kinetic and internal energy for compressible runs,
! including the mass density
!
! Output file contains:
! 'compr_energy.txt': time, kinetic energy, internal energy
!
! Parameters
!     gam1: gamma-1 constant (adiabatic constant)
!     cp1 : 2/(gam1 smach^2), the coefficient of the enthalpy
!     d   : input matrix with density (in Fourier space)
!     a   : input matrix with v_x (in Fourier space)
!     b   : input matrix with v_y (in Fourier space)
!     c   : input matrix with v_z (in Fourier space) [A = (a,b,c)]
!     t   : number of time steps made
!     dt  : time step
!     path: path for the output
!
      USE fprecision
      USE kes
      USE grid
      USE commtypes
      USE mpivars
      USE fft
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: b
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: d
      REAL(KIND=GP),    INTENT(IN)  :: gam1,cp1,dt
      INTEGER,          INTENT(IN)  :: t
      CHARACTER(len=*), INTENT(IN)  :: path
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: e
      REAL(KIND=GP)   , POINTER, DIMENSION(:,:,:) :: r1,r4
      REAL(KIND=GP)                 :: tmp, tmp1, gam0
      DOUBLE PRECISION              :: tot_ekin
      DOUBLE PRECISION              :: loc_ekin
      DOUBLE PRECISION              :: tot_eint
      DOUBLE PRECISION              :: loc_eint
      INTEGER                       :: i,j,k,m
      LOGICAL                       :: bret

      CALL gws%get_complex_htmp(e,bret)
      CALL gws%get_real_htmp(r1,bret)
      CALL gws%get_real_htmp(r4,bret)

      tot_ekin = 0.0D0
      tot_eint = 0.0D0
      loc_ekin = 0.0D0
      loc_eint = 0.0D0
      tmp  = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**3
      tmp1 = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
      gam0 = gam1 + 1.0_GP

      CALL copy3(d,e)
      CALL fftp3d_complex_to_real(plancr,e,r4,MPI_COMM_WORLD)
!$omp parallel do collapse(2) private (i) reduction(+:loc_eint)
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               loc_eint = loc_eint + (r4(i,j,k)*tmp1)**gam0
            END DO
         END DO
      END DO

      ! Kinetic energy: one velocity component at a time
      DO m = 1,3
         IF (m.eq.1) THEN
            CALL copy3(a,e)
         ELSE IF (m.eq.2) THEN
            CALL copy3(b,e)
         ELSE
            CALL copy3(c,e)
         ENDIF
         CALL fftp3d_complex_to_real(plancr,e,r1,MPI_COMM_WORLD)
!$omp parallel do collapse(2) private (i) reduction(+:loc_ekin)
         DO k = ksta,kend
            DO j = 1,ny
               DO i = 1,nx
                  loc_ekin = loc_ekin + r4(i,j,k)*r1(i,j,k)*r1(i,j,k)*tmp
               END DO
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
         OPEN(1,file=trim(path) // '/compr_energy.txt',position='append')
         WRITE(1,FMT='(E13.6,E22.14,E22.14)') (t-1)*dt,tot_ekin,tot_eint
         CLOSE(1)
      ENDIF

      CALL gws%free_real_htmp(r4)
      CALL gws%free_real_htmp(r1)
      CALL gws%free_complex_htmp(e)

      RETURN
      END SUBROUTINE energycompr

!*****************************************************************
      SUBROUTINE massenergycompi(gam1,a,b,c,d,e,t,dt,path)
!-----------------------------------------------------------------
!
! Computes and outputs the total mass and the kinetic and
! internal energy densities for compressible runs, for the
! system of PDEs in which the internal energy is evolved.
!
! Output file contains:
! 'compi_massenergy.txt': time, mass, kinetic energy density,
! internal energy density, rms velocity, rms sound speed, mean
! Mach number
!
! Parameters
!     gam1: gamma - 1
!     a   : input matrix with v_x (in Fourier space)
!     b   : input matrix with v_y (in Fourier space)
!     c   : input matrix with v_z (in Fourier space) [A = (a,b,c)]
!     d   : input matrix with density (in Fourier space)
!     e   : input matrix with internal energy density (in Fourier space)
!     t   : number of time steps made
!     dt  : time step
!     path: path for the output
!
      USE fprecision
      USE kes
      USE grid
      USE commtypes
      USE mpivars
      USE fft
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: b
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: d
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: e
      REAL(KIND=GP),    INTENT(IN)  :: gam1, dt
      INTEGER,          INTENT(IN)  :: t
      CHARACTER(len=*), INTENT(IN)  :: path
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: t1
      REAL(KIND=GP)   , POINTER, DIMENSION(:,:,:) :: r1,r4,r5
      REAL(KIND=GP)                 :: csq, tmp1, tmp2, tmp3, vsq
      DOUBLE PRECISION              :: tot_ekin,tot_eint,tot_mass,tot_mach
      DOUBLE PRECISION              :: tot_c,tot_v
      DOUBLE PRECISION              :: tiny,v2,vloc(5),vtot(5)
      INTEGER                       :: i,j,k,m
      LOGICAL                       :: bret

      CALL gws%get_complex_htmp(t1,bret)
      CALL gws%get_real_htmp(r1,bret)
      CALL gws%get_real_htmp(r4,bret)
      CALL gws%get_real_htmp(r5,bret)

      tiny = 100.0*epsilon(tot_mach)
      vloc = 0.0D0
      tmp3  = 1.0_GP/ &
              (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**3
      tmp2  = 1.0_GP/ &
              (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      tmp1 = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))

      ! v^2 in r5, accumulated one component at a time
      DO m = 1,3
         IF (m.eq.1) THEN
            CALL copy3(a,t1)
         ELSE IF (m.eq.2) THEN
            CALL copy3(b,t1)
         ELSE
            CALL copy3(c,t1)
         ENDIF
         CALL fftp3d_complex_to_real(plancr,t1,r1,MPI_COMM_WORLD)
         IF (m.eq.1) THEN
!$omp parallel do collapse(2) private (i)
            DO k = ksta,kend
               DO j = 1,ny
                  DO i = 1,nx
                     r5(i,j,k) = r1(i,j,k)*r1(i,j,k)
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do collapse(2) private (i)
            DO k = ksta,kend
               DO j = 1,ny
                  DO i = 1,nx
                     r5(i,j,k) = r5(i,j,k) + r1(i,j,k)*r1(i,j,k)
                  END DO
               END DO
            END DO
         ENDIF
      END DO
      ! Density in r4, internal energy in r1
      CALL copy3(d,t1)
      CALL fftp3d_complex_to_real(plancr,t1,r4,MPI_COMM_WORLD)
      CALL copy3(e,t1)
      CALL fftp3d_complex_to_real(plancr,t1,r1,MPI_COMM_WORLD)

!$omp parallel do collapse(2) private (i,v2,vsq,csq) reduction(+:vloc)
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               v2      = r5(i,j,k)
               vsq     = v2 * tmp2
               csq     = gam1*(gam1+1.0_GP) * r1(i,j,k) / (r4(i,j,k)+tiny)
               vloc(1) = vloc(1) + (r4(i,j,k) * v2 * tmp3)
               vloc(2) = vloc(2) + (r1(i,j,k)*tmp1)
               vloc(3) = vloc(3) + (r4(i,j,k)*tmp1)
               vloc(4) = vloc(4) + vsq
               vloc(5) = vloc(5) + csq
            END DO
         END DO
      END DO

      ! Compute averages over grid:
      vloc = vloc*tmp1

      CALL MPI_REDUCE(vloc,vtot,5,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      tot_ekin = vtot(1)
      tot_eint = vtot(2)
      tot_mass = vtot(3)
      tot_v    = sqrt(vtot(4))
      tot_c    = sqrt(vtot(5))
      tot_mach = tot_v/(tot_c+tiny)

      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(path) // '/compi_massenergy.txt',position='append')
         WRITE(1,10) (t-1)*dt,tot_mass,tot_ekin,tot_eint,tot_v,tot_c,tot_mach
   10    FORMAT( E13.6,E26.18,E26.18,E26.18,E26.18,E26.18,E26.18 )
         CLOSE(1)
      ENDIF

      CALL gws%free_real_htmp(r5)
      CALL gws%free_real_htmp(r4)
      CALL gws%free_real_htmp(r1)
      CALL gws%free_complex_htmp(t1)

      RETURN
      END SUBROUTINE massenergycompi

END MODULE pseudospec_compr
