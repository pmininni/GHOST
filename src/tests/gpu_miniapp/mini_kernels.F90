!=================================================================
! mini_kernels.F90
! Field kernels of the mini-app in the dual-pragma form (OpenMP
! threads on the host, OpenMP target offload on the GPU), plus
! serial host reference versions (_ref) used to check them.
!
! The loop headers are written out in full in every kernel rather
! than pulled from an include file: flang identifies offload entries
! by source location, and several target regions expanded from the
! same include line inside one routine get their device kernels
! swapped (observed with ROCm 7.2 amdflang in derivk3).
!
! The kernels receive explicit-shape dummy arrays, as the GHOST
! pseudospectral routines do. In an offload build the arrays are
! expected to be resident on the device already (mapped at their
! allocation point); the implicit map of the dummies then finds
! them in the device table and no data is transferred.
!=================================================================
#include "gpu_defs.h"

  MODULE mini_kernels
      USE mini_prec
      USE mini_grid
      USE mini_mpivars
      USE mini_kes
      USE mini_state
      IMPLICIT NONE
      COMPLEX(KIND=GP), PARAMETER :: im = (0.0_GP,1.0_GP)

    CONTAINS

!*****************************************************************
      SUBROUTINE derivk3(a,b,dir)
!
! Derivative in Fourier space, dir = 1,2,3 (as in pseudospec3D_hd)
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: b
      INTEGER, INTENT(IN) :: dir
      INTEGER             :: i,j,k

      IF (dir.eq.1) THEN
#if defined(GHOST_GPU)
#if GPU_LOOP_SPELLING == 2
!$omp target teams loop collapse(3)
#else
!$omp target teams distribute parallel do collapse(3)
#endif
      DO i = ista,iend
        DO j = 1,ny
          DO k = 1,nz
#else
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
        DO j = 1,ny
          DO CONCURRENT (k=1:nz)
#endif
                  b(k,j,i) = im*kx(i)*a(k,j,i)
               END DO
            END DO
         END DO
      ELSE IF (dir.eq.2) THEN
#if defined(GHOST_GPU)
#if GPU_LOOP_SPELLING == 2
!$omp target teams loop collapse(3)
#else
!$omp target teams distribute parallel do collapse(3)
#endif
      DO i = ista,iend
        DO j = 1,ny
          DO k = 1,nz
#else
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
        DO j = 1,ny
          DO CONCURRENT (k=1:nz)
#endif
                  b(k,j,i) = im*ky(j)*a(k,j,i)
               END DO
            END DO
         END DO
      ELSE
#if defined(GHOST_GPU)
#if GPU_LOOP_SPELLING == 2
!$omp target teams loop collapse(3)
#else
!$omp target teams distribute parallel do collapse(3)
#endif
      DO i = ista,iend
        DO j = 1,ny
          DO k = 1,nz
#else
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
        DO j = 1,ny
          DO CONCURRENT (k=1:nz)
#endif
                  b(k,j,i) = im*kz(k)*a(k,j,i)
               END DO
            END DO
         END DO
      ENDIF
      END SUBROUTINE derivk3

!*****************************************************************
      SUBROUTINE laplak3(a,b)
!
! Laplacian in Fourier space, uses the module array kk2
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: b
      INTEGER :: i,j,k

#if defined(GHOST_GPU)
#if GPU_LOOP_SPELLING == 2
!$omp target teams loop collapse(3)
#else
!$omp target teams distribute parallel do collapse(3)
#endif
      DO i = ista,iend
        DO j = 1,ny
          DO k = 1,nz
#else
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
        DO j = 1,ny
          DO CONCURRENT (k=1:nz)
#endif
               b(k,j,i) = -kk2(k,j,i)*a(k,j,i)
            END DO
         END DO
      END DO
      END SUBROUTINE laplak3

!*****************************************************************
      SUBROUTINE axpy3(y,x,a)
!
! y = x + a*y on explicit-shape arrays (variant "argument" of the
! stepper update: the state components are passed as arguments)
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: y
      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(nz,ny,ista:iend) :: x
      REAL(KIND=GP), INTENT(IN) :: a
      INTEGER :: i,j,k

#if defined(GHOST_GPU)
#if GPU_LOOP_SPELLING == 2
!$omp target teams loop collapse(3)
#else
!$omp target teams distribute parallel do collapse(3)
#endif
      DO i = ista,iend
        DO j = 1,ny
          DO k = 1,nz
#else
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
        DO j = 1,ny
          DO CONCURRENT (k=1:nz)
#endif
               y(k,j,i) = x(k,j,i) + a*y(k,j,i)
            END DO
         END DO
      END DO
      END SUBROUTINE axpy3

!*****************************************************************
      SUBROUTINE update_arg(uout,uin,nstate,dt)
!
! Stepper update, variant "argument": one kernel call per component
      TYPE(GStateComp), INTENT(INOUT) :: uout(:)
      TYPE(GStateComp), INTENT   (IN) :: uin(:)
      INTEGER, INTENT(IN) :: nstate
      REAL(KIND=GP), INTENT(IN) :: dt
      INTEGER :: ic

      DO ic = 1,nstate
         CALL axpy3(uout(ic)%ccomp,uin(ic)%ccomp,dt)
      END DO
      END SUBROUTINE update_arg

!*****************************************************************
      SUBROUTINE update_ptr(uout,uin,nstate,dt)
!
! Stepper update, variant "pointer": the kernel references local
! pointers aliased to the state components
      TYPE(GStateComp), INTENT(INOUT), TARGET :: uout(:)
      TYPE(GStateComp), INTENT   (IN), TARGET :: uin(:)
      INTEGER, INTENT(IN) :: nstate
      REAL(KIND=GP), INTENT(IN) :: dt
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: po,pi
      INTEGER :: ic,i,j,k

      DO ic = 1,nstate
         po => uout(ic)%ccomp
         pi => uin (ic)%ccomp
#if defined(GHOST_GPU)
#if GPU_LOOP_SPELLING == 2
!$omp target teams loop collapse(3)
#else
!$omp target teams distribute parallel do collapse(3)
#endif
      DO i = ista,iend
        DO j = 1,ny
          DO k = 1,nz
#else
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
        DO j = 1,ny
          DO CONCURRENT (k=1:nz)
#endif
                  po(k,j,i) = pi(k,j,i) + dt*po(k,j,i)
               END DO
            END DO
         END DO
      END DO
      END SUBROUTINE update_ptr

#if defined(TEST_DIRECT_DT)
!*****************************************************************
      SUBROUTINE update_direct(uout,uin,nstate,dt)
!
! Stepper update, variant "direct": the kernel indexes the derived
! type components directly, as canuto_stepper.f90 does today
      TYPE(GStateComp), INTENT(INOUT) :: uout(:)
      TYPE(GStateComp), INTENT   (IN) :: uin(:)
      INTEGER, INTENT(IN) :: nstate
      REAL(KIND=GP), INTENT(IN) :: dt
      INTEGER :: ic,i,j,k

#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(4)
      DO ic = 1,nstate
        DO i = ista,iend
          DO j = 1,ny
            DO k = 1,nz
#else
!$omp parallel do collapse(3) private (k)
      DO ic = 1,nstate
        DO i = ista,iend
          DO j = 1,ny
            DO CONCURRENT (k=1:nz)
#endif
              uout(ic)%ccomp(k,j,i) = uin(ic)%ccomp(k,j,i) + dt*uout(ic)%ccomp(k,j,i)
            END DO
          END DO
        END DO
      END DO
      END SUBROUTINE update_direct
#endif

!*****************************************************************
      SUBROUTINE energy_dev(a,b,c,ek)
!
! Local contribution to the energy, with a reduction inside the
! kernel. Hermitian weight 2 for kx > 0 as in GHOST's energy. The
! squared modulus is written out as REAL**2+AIMAG**2: ABS of a
! complex inside a kernel becomes a call with a complex argument,
! which upstream flang (23.1) cannot yet lower for NVPTX.
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      REAL(KIND=GP), INTENT(OUT) :: ek
      REAL(KIND=GP) :: tmp
      INTEGER :: i,j,k

      tmp = 0.0_GP
#if defined(GHOST_GPU)
#if GPU_LOOP_SPELLING == 2
!$omp target teams loop collapse(3) reduction(+:tmp)
#else
!$omp target teams distribute parallel do collapse(3) reduction(+:tmp)
#endif
#else
!$omp parallel do collapse(3) reduction(+:tmp)
#endif
      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               IF (i.eq.1) THEN
                  tmp = tmp + (REAL(a(k,j,i))**2+AIMAG(a(k,j,i))**2 + &
                               REAL(b(k,j,i))**2+AIMAG(b(k,j,i))**2 + &
                               REAL(c(k,j,i))**2+AIMAG(c(k,j,i))**2)
               ELSE
                  tmp = tmp + 2*(REAL(a(k,j,i))**2+AIMAG(a(k,j,i))**2 + &
                                 REAL(b(k,j,i))**2+AIMAG(b(k,j,i))**2 + &
                                 REAL(c(k,j,i))**2+AIMAG(c(k,j,i))**2)
               ENDIF
            END DO
         END DO
      END DO
      ek = tmp
      END SUBROUTINE energy_dev

!*****************************************************************
      SUBROUTINE rhs_mask(lap,nl,f,out,nu)
!
! Dealiased right-hand side assembly, as the last loop of dudt
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: lap,nl,f
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: out
      REAL(KIND=GP), INTENT(IN) :: nu
      INTEGER :: i,j,k

#if defined(GHOST_GPU)
#if GPU_LOOP_SPELLING == 2
!$omp target teams loop collapse(3)
#else
!$omp target teams distribute parallel do collapse(3)
#endif
      DO i = ista,iend
        DO j = 1,ny
          DO k = 1,nz
#else
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
        DO j = 1,ny
          DO CONCURRENT (k=1:nz)
#endif
               IF ((kk2(k,j,i).le.kmax).and.(kk2(k,j,i).ge.tiny)) THEN
                  out(k,j,i) = nu*lap(k,j,i) + nl(k,j,i) + f(k,j,i)
               ELSE
                  out(k,j,i) = 0.0_GP
               ENDIF
            END DO
         END DO
      END DO
      END SUBROUTINE rhs_mask

!*****************************************************************
      SUBROUTINE rprod(a,b,c)
!
! Pointwise product in real space (a nonlinear term)
      REAL(KIND=GP), INTENT (IN), DIMENSION(nx,ny,ksta:kend) :: a,b
      REAL(KIND=GP), INTENT(OUT), DIMENSION(nx,ny,ksta:kend) :: c
      INTEGER :: i,j,k

#if defined(GHOST_GPU)
#if GPU_LOOP_SPELLING == 2
!$omp target teams loop collapse(3)
#else
!$omp target teams distribute parallel do collapse(3)
#endif
      DO k = ksta,kend
        DO j = 1,ny
          DO i = 1,nx
#else
!$omp parallel do collapse(2) private (i)
      DO k = ksta,kend
        DO j = 1,ny
          DO CONCURRENT (i=1:nx)
#endif
               c(i,j,k) = a(i,j,k)*b(i,j,k)
            END DO
         END DO
      END DO
      END SUBROUTINE rprod

!*****************************************************************
      SUBROUTINE probe_k(b)
!
! Diagnostic: writes the device copies of the wavenumber arrays
! into b, b(k,j,i) = (kx(i)+ky(j)) + i kz(k)
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: b
      INTEGER :: i,j,k

#if defined(GHOST_GPU)
#if GPU_LOOP_SPELLING == 2
!$omp target teams loop collapse(3)
#else
!$omp target teams distribute parallel do collapse(3)
#endif
      DO i = ista,iend
        DO j = 1,ny
          DO k = 1,nz
#else
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
        DO j = 1,ny
          DO CONCURRENT (k=1:nz)
#endif
               b(k,j,i) = CMPLX(kx(i)+ky(j),kz(k),KIND=GP)
            END DO
         END DO
      END DO
      END SUBROUTINE probe_k

!*****************************************************************
      SUBROUTINE probe_im(a,b)
!
! Diagnostic: complex multiplication by the imaginary unit
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: b
      INTEGER :: i,j,k

#if defined(GHOST_GPU)
#if GPU_LOOP_SPELLING == 2
!$omp target teams loop collapse(3)
#else
!$omp target teams distribute parallel do collapse(3)
#endif
      DO i = ista,iend
        DO j = 1,ny
          DO k = 1,nz
#else
!$omp parallel do collapse(2) private (k)
      DO i = ista,iend
        DO j = 1,ny
          DO CONCURRENT (k=1:nz)
#endif
               b(k,j,i) = im*a(k,j,i)
            END DO
         END DO
      END DO
      END SUBROUTINE probe_im

! ================================================================
! Serial host reference versions
! ================================================================
      SUBROUTINE derivk3_ref(a,b,dir)
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: b
      INTEGER, INTENT(IN) :: dir
      INTEGER :: i,j,k
      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               IF (dir.eq.1) THEN
                  b(k,j,i) = im*kx(i)*a(k,j,i)
               ELSE IF (dir.eq.2) THEN
                  b(k,j,i) = im*ky(j)*a(k,j,i)
               ELSE
                  b(k,j,i) = im*kz(k)*a(k,j,i)
               ENDIF
            END DO
         END DO
      END DO
      END SUBROUTINE derivk3_ref

      SUBROUTINE laplak3_ref(a,b)
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: b
      INTEGER :: i,j,k
      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               b(k,j,i) = -kk2(k,j,i)*a(k,j,i)
            END DO
         END DO
      END DO
      END SUBROUTINE laplak3_ref

      SUBROUTINE energy_ref(a,b,c,ek)
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      REAL(KIND=GP), INTENT(OUT) :: ek
      DOUBLE PRECISION :: tmp
      INTEGER :: i,j,k
      tmp = 0.0D0
      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               IF (i.eq.1) THEN
                  tmp = tmp + (REAL(a(k,j,i))**2+AIMAG(a(k,j,i))**2 + &
                               REAL(b(k,j,i))**2+AIMAG(b(k,j,i))**2 + &
                               REAL(c(k,j,i))**2+AIMAG(c(k,j,i))**2)
               ELSE
                  tmp = tmp + 2*(REAL(a(k,j,i))**2+AIMAG(a(k,j,i))**2 + &
                                 REAL(b(k,j,i))**2+AIMAG(b(k,j,i))**2 + &
                                 REAL(c(k,j,i))**2+AIMAG(c(k,j,i))**2)
               ENDIF
            END DO
         END DO
      END DO
      ek = real(tmp,kind=GP)
      END SUBROUTINE energy_ref

      SUBROUTINE rhs_mask_ref(lap,nl,f,out,nu)
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: lap,nl,f
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: out
      REAL(KIND=GP), INTENT(IN) :: nu
      INTEGER :: i,j,k
      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               IF ((kk2(k,j,i).le.kmax).and.(kk2(k,j,i).ge.tiny)) THEN
                  out(k,j,i) = nu*lap(k,j,i) + nl(k,j,i) + f(k,j,i)
               ELSE
                  out(k,j,i) = 0.0_GP
               ENDIF
            END DO
         END DO
      END DO
      END SUBROUTINE rhs_mask_ref

      SUBROUTINE rprod_ref(a,b,c)
      REAL(KIND=GP), INTENT (IN), DIMENSION(nx,ny,ksta:kend) :: a,b
      REAL(KIND=GP), INTENT(OUT), DIMENSION(nx,ny,ksta:kend) :: c
      INTEGER :: i,j,k
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               c(i,j,k) = a(i,j,k)*b(i,j,k)
            END DO
         END DO
      END DO
      END SUBROUTINE rprod_ref

  END MODULE mini_kernels
