!=================================================================
! mini_main.F90
! Driver of the GHOST GPU mini-app. Validates, per compiler and
! vendor, the assumptions of the offload design:
!   1. dual-pragma kernels on explicit-shape dummies that are
!      resident on the device (no transfers inside the step),
!   2. the three ways of updating an array of derived types
!      (direct component indexing, pointer alias, argument),
!   3. workspace pool pointers resolving in the device table,
!   4. device reductions and their run-to-run determinism,
!   5. cuFFT/hipFFT on resident arrays through use_device_addr,
!   6. the slab exchange with GPU-aware MPI (datatypes, packed) and
!      with host staging, with timings.
! Every device result is compared against a host reference.
!
! Usage: gpu_miniapp [N] [niter]
!=================================================================
#include "gpu_defs.h"

      PROGRAM gpu_miniapp
      USE ISO_C_BINDING
      USE mini_prec
      USE mini_mpivars
      USE mini_grid
      USE mini_kes
      USE mini_state
      USE mini_workspace
      USE mini_device
      USE mini_xfer
      USE mini_kernels
      USE mini_fft
      IMPLICIT NONE

      TYPE(GStateComp), ALLOCATABLE, TARGET :: field(:),fnxt(:)
      TYPE(GWorkspace) :: ws
      TYPE(FFTPLAN)    :: planrc,plancr
      COMPLEX(KIND=GP), POINTER :: C1(:,:,:),C2(:,:,:),C3(:,:,:)
      REAL(KIND=GP)   , POINTER :: R1(:,:,:),R2(:,:,:),R3(:,:,:)
      COMPLEX(KIND=GP), ALLOCATABLE :: href(:,:,:)
      REAL(KIND=GP)   , ALLOCATABLE :: rref(:,:,:)
      REAL(KIND=GP)    :: dt,nu,ek1,ek2,ekr,x,y,z,tmp
      DOUBLE PRECISION :: t0,t1,tt(4),ttmax(4),tt3(3),tt3max(3)
      INTEGER :: n,niter,ic,dir,mode,it,i,j,k,nfail
      LOGICAL :: bret,ondev,ok
      CHARACTER(len=32) :: arg

! -------------------- Initialization ----------------------------
      CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      n = 64; niter = 10
      IF (COMMAND_ARGUMENT_COUNT().ge.1) THEN
         CALL GET_COMMAND_ARGUMENT(1,arg); READ(arg,*) n
      ENDIF
      IF (COMMAND_ARGUMENT_COUNT().ge.2) THEN
         CALL GET_COMMAND_ARGUMENT(2,arg); READ(arg,*) niter
      ENDIF
      nx = n; ny = n; nz = n
      CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
      CALL range(1,nz    ,nprocs,myrank,ksta,kend)
#if defined(GDOUBLE_PRECISION)
      GC_REAL = MPI_DOUBLE_PRECISION; GC_COMPLEX = MPI_DOUBLE_COMPLEX
#else
      GC_REAL = MPI_REAL; GC_COMPLEX = MPI_COMPLEX
#endif
      CALL device_init(myrank)
#if defined(GHOST_GPU)
      ondev = .TRUE.
#else
      ondev = .FALSE.
#endif
      nfail = 0

      IF (myrank.eq.0) THEN
         PRINT '(A)'      ,'GHOST GPU mini-app'
         PRINT '(A,I6)'   ,'  grid N        = ',n
         PRINT '(A,I6)'   ,'  MPI ranks     = ',nprocs
         PRINT '(A,I6)'   ,'  devices/node  = ',ndev
#if defined(GDOUBLE_PRECISION)
         PRINT '(A)'      ,'  precision     = DOUBLE'
#else
         PRINT '(A)'      ,'  precision     = SINGLE'
#endif
#if defined(GPU_NVIDIA)
         PRINT '(A)'      ,'  build         = offload, NVIDIA (cuFFT)'
#elif defined(GPU_AMD)
         PRINT '(A)'      ,'  build         = offload, AMD (hipFFT)'
#else
         PRINT '(A)'      ,'  build         = host OpenMP (FFTW)'
#endif
#if defined(GHOST_GPU)
         PRINT '(A,I2)'   ,'  loop spelling = ',GPU_LOOP_SPELLING
#endif
      ENDIF
      PRINT '(A,I3,A,I3,A,I5,A,I5,A,I5,A,I5)','  rank ',myrank,' device ',&
            targetdev,'  kx ',ista,':',iend,'  z ',ksta,':',kend

      CALL kes_init()
      CALL GState_alloc(field,3)
      CALL GState_alloc(fnxt ,3)
      CALL ws%initialize_pool(3,4)
      CALL fftp3d_create_plan(planrc,(/nx,ny,nz/),FFT_FORWARD)
      CALL fftp3d_create_plan(plancr,(/nx,ny,nz/),FFT_BACKWARD)
      ALLOCATE( href(nz,ny,ista:iend), rref(nx,ny,ksta:kend) )

! -------------------- Initial condition -------------------------
! Taylor-Green flow in real space, transformed with the host path,
! then copied to the device. field(3) stays zero.
      CALL ws%get_real_tmp(R1,bret)
      CALL fill_tg(R1,1)
      CALL fftp3d_real_to_complex(planrc,R1,field(1)%ccomp,MPI_COMM_WORLD,.FALSE.)
      CALL fill_tg(R1,2)
      CALL fftp3d_real_to_complex(planrc,R1,field(2)%ccomp,MPI_COMM_WORLD,.FALSE.)
      CALL ws%free_real_tmp(R1)
      CALL GState_update_to(field)
      DO ic = 1,3
         fnxt(ic)%ccomp = field(ic)%ccomp
      END DO
      CALL GState_update_to(fnxt)
      dt = 0.01_GP
      nu = 0.1_GP

      IF (myrank.eq.0) PRINT '(/A)','--- Kernels on resident arrays ---'

! -------------------- T1: derivk3 -------------------------------
      CALL ws%get_complex_tmp(C1,bret)
      DO dir = 1,3
         CALL derivk3(field(1)%ccomp,C1,dir)
         CALL cupdate_from(C1)
         CALL derivk3_ref(field(1)%ccomp,href,dir)
         CALL check_c('derivk3 dir='//ACHAR(48+dir),C1,href,1e-5_GP)
      END DO

! -------------------- T1b: diagnostics for derivk3 --------------
      CALL probe_k(C1)
      CALL cupdate_from(C1)
      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               href(k,j,i) = CMPLX(kx(i)+ky(j),kz(k),KIND=GP)
            END DO
         END DO
      END DO
      CALL check_c('probe: device kx,ky,kz',C1,href,1e-6_GP)
      CALL probe_im(field(1)%ccomp,C1)
      CALL cupdate_from(C1)
      href = im*field(1)%ccomp
      CALL check_c('probe: im*a',C1,href,1e-6_GP)

! -------------------- T2: laplak3 (module array kk2) ------------
      CALL laplak3(field(2)%ccomp,C1)
      CALL cupdate_from(C1)
      CALL laplak3_ref(field(2)%ccomp,href)
      CALL check_c('laplak3',C1,href,1e-5_GP)

! -------------------- T3: rhs with dealiasing mask --------------
      CALL ws%get_complex_tmp(C2,bret)
      CALL ws%get_complex_tmp(C3,bret)
      CALL laplak3(field(1)%ccomp,C2)
      CALL derivk3(field(2)%ccomp,C3,3)
      CALL rhs_mask(C2,C3,field(3)%ccomp,C1,nu)
      CALL cupdate_from(C1); CALL cupdate_from(C2); CALL cupdate_from(C3)
      CALL rhs_mask_ref(C2,C3,field(3)%ccomp,href,nu)
      CALL check_c('rhs_mask',C1,href,1e-5_GP)

! -------------------- T4: pool pointers present on device -------
      ok = device_present(C_LOC(C1)) .and. device_present(C_LOC(field(1)%ccomp))
      IF (myrank.eq.0) THEN
         IF (ok) THEN
            PRINT '(A,T44,A)','present(pool ptr, state comp)','PASS'
         ELSE
            PRINT '(A,T44,A)','present(pool ptr, state comp)','FAIL'
            nfail = nfail+1
         ENDIF
      ENDIF

! -------------------- T5: real-space product --------------------
      CALL ws%get_real_tmp(R1,bret)
      CALL ws%get_real_tmp(R2,bret)
      CALL ws%get_real_tmp(R3,bret)
      CALL fill_tg(R1,1); CALL fill_tg(R2,2)
      CALL rupdate_to(R1); CALL rupdate_to(R2)
      CALL rprod(R1,R2,R3)
      CALL rupdate_from(R3)
      CALL rprod_ref(R1,R2,rref)
      CALL check_r('rprod (real space)',R3,rref,1e-6_GP)
      CALL ws%free_real_tmp(R3)

      IF (myrank.eq.0) PRINT '(/A)','--- Stepper update on arrays of derived types ---'

! -------------------- T6: update variants -----------------------
! fnxt = field + dt*fnxt with fnxt = field at input: (1+dt)*field
      href = (1.0_GP+dt)*field(1)%ccomp
      CALL update_arg(fnxt,field,3,dt)
      CALL cupdate_from(fnxt(1)%ccomp)
      CALL check_c('update, argument variant',fnxt(1)%ccomp,href,1e-6_GP)
      CALL reset_fnxt()

      CALL update_ptr(fnxt,field,3,dt)
      CALL cupdate_from(fnxt(1)%ccomp)
      CALL check_c('update, pointer variant',fnxt(1)%ccomp,href,1e-6_GP)
      CALL reset_fnxt()

#if defined(TEST_DIRECT_DT)
      CALL update_direct(fnxt,field,3,dt)
      CALL cupdate_from(fnxt(1)%ccomp)
      CALL check_c('update, direct component variant',fnxt(1)%ccomp,href,1e-6_GP)
      CALL reset_fnxt()
#else
      IF (myrank.eq.0) PRINT '(A,T44,A)','update, direct component variant','NOT COMPILED'
#endif

      IF (myrank.eq.0) PRINT '(/A)','--- Reductions ---'

! -------------------- T7: energy reduction ----------------------
      CALL energy_dev(field(1)%ccomp,field(2)%ccomp,field(3)%ccomp,ek1)
      CALL energy_dev(field(1)%ccomp,field(2)%ccomp,field(3)%ccomp,ek2)
      CALL energy_ref(field(1)%ccomp,field(2)%ccomp,field(3)%ccomp,ekr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,ek1,1,GC_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,ek2,1,GC_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,ekr,1,GC_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         tmp = abs(ek1-ekr)/abs(ekr)
         CALL verdict('energy reduction vs host',tmp,1e-5_GP)
         IF (ek1.eq.ek2) THEN
            PRINT '(A,T44,A)','energy reduction run-to-run','BITWISE EQUAL'
         ELSE
            PRINT '(A,T44,A,ES10.2)','energy reduction run-to-run','DIFFERS, rel ', &
                  abs(ek1-ek2)/abs(ekr)
         ENDIF
      ENDIF

      IF (myrank.eq.0) PRINT '(/A)','--- Parallel FFT: device vs host (FFTW) ---'

! -------------------- T8: forward transform, each exchange mode -
! R1 holds the TG field on host and device. href is the host result.
      CALL fftp3d_real_to_complex(planrc,R1,href,MPI_COMM_WORLD,.FALSE.)
      DO mode = 1,3
         IF (.not.ondev .and. mode.eq.3) EXIT
         exchange_mode = mode
         CALL fftp3d_real_to_complex(planrc,R1,C1,MPI_COMM_WORLD,ondev)
         CALL cupdate_from(C1)
         CALL check_c('r2c, exchange mode '//ACHAR(48+mode),C1,href,1e-4_GP)
         CALL fftp3d_complex_to_real(plancr,C1,R2,MPI_COMM_WORLD,ondev)
         CALL rupdate_from(R2)
         tmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
         rref = R2*tmp
         CALL check_r('c2r roundtrip, mode '//ACHAR(48+mode),rref,R1,1e-4_GP)
      END DO

! -------------------- T9: timings -------------------------------
      IF (myrank.eq.0) THEN
         PRINT '(/A,I4,A)','--- Timings per r2c+c2r pair, max over ranks, ',niter,' iterations ---'
         PRINT '(A)','path                     total[ms]    fft[ms]   comm[ms]  transp[ms]'
      ENDIF
      DO mode = 0,3
         IF (.not.ondev .and. mode.eq.3) EXIT
         IF (mode.eq.0) THEN
            exchange_mode = 1
         ELSE
            exchange_mode = mode
         ENDIF
         CALL fft_timers_reset()
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         t0 = MPI_WTIME()
         DO it = 1,niter
            CALL fftp3d_real_to_complex(planrc,R1,C1,MPI_COMM_WORLD,ondev.and.(mode.gt.0))
            CALL fftp3d_complex_to_real(plancr,C1,R2,MPI_COMM_WORLD,ondev.and.(mode.gt.0))
         END DO
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         t1 = MPI_WTIME()
         tt = (/ (t1-t0), tfft, tcom, ttra /)*1000.0D0/niter
         CALL MPI_REDUCE(tt,ttmax,4,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
         tt3 = (/ tpack, tself, tmpi /)*1000.0D0/niter
         CALL MPI_REDUCE(tt3,tt3max,3,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) THEN
            IF (mode.eq.0) THEN
               PRINT '(A,T24,4F11.3)','host (FFTW, types)',ttmax
            ELSE IF (mode.eq.1) THEN
               PRINT '(A,T24,4F11.3)','device, types (1)',ttmax
            ELSE IF (mode.eq.2) THEN
               PRINT '(A,T24,4F11.3)','device, packed (2)',ttmax
               PRINT '(A,T24,3F11.3)','   (2) pack/self/mpi',tt3max
            ELSE
               PRINT '(A,T24,4F11.3)','device, host staged (3)',ttmax
            ENDIF
         ENDIF
      END DO

! -------------------- T10: kernel timings -----------------------
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t0 = MPI_WTIME()
      DO it = 1,niter
         CALL derivk3(field(1)%ccomp,C1,1)
         CALL laplak3(field(2)%ccomp,C2)
         CALL rhs_mask(C2,C1,field(3)%ccomp,C3,nu)
      END DO
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      t1 = MPI_WTIME()
      tt(1) = (t1-t0)*1000.0D0/niter
      CALL MPI_REDUCE(tt,ttmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) PRINT '(A,T24,F11.3)','3 kernels (deriv,lap,rhs)',ttmax(1)

! -------------------- Summary -----------------------------------
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,nfail,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF (nfail.eq.0) THEN
            PRINT '(/A)','ALL TESTS PASSED'
         ELSE
            PRINT '(/A,I3,A)','FAILED: ',nfail,' checks'
         ENDIF
      ENDIF

      CALL ws%free_real_tmp(R1)
      CALL ws%free_real_tmp(R2)
      CALL ws%free_complex_tmp(C1)
      CALL ws%free_complex_tmp(C2)
      CALL ws%free_complex_tmp(C3)
      CALL fftp3d_destroy_plan(planrc)
      CALL fftp3d_destroy_plan(plancr)
      CALL ws%destroy_pool()
      CALL GState_dealloc(field)
      CALL GState_dealloc(fnxt)
      DEALLOCATE( href,rref )
      CALL MPI_FINALIZE(ierr)

    CONTAINS

!-----------------------------------------------------------------
      SUBROUTINE fill_tg(r,comp)
!
! Taylor-Green velocity component in real space (host only)
      REAL(KIND=GP), INTENT(OUT) :: r(nx,ny,ksta:kend)
      INTEGER, INTENT(IN) :: comp
      REAL(KIND=GP), PARAMETER :: twopi = 6.283185307179586_GP
      INTEGER :: i,j,k
      DO k = ksta,kend
         z = twopi*real(k-1,kind=GP)/real(nz,kind=GP)
         DO j = 1,ny
            y = twopi*real(j-1,kind=GP)/real(ny,kind=GP)
            DO i = 1,nx
               x = twopi*real(i-1,kind=GP)/real(nx,kind=GP)
               IF (comp.eq.1) THEN
                  r(i,j,k) =  sin(x)*cos(y)*cos(z)
               ELSE
                  r(i,j,k) = -cos(x)*sin(y)*cos(z)
               ENDIF
            END DO
         END DO
      END DO
      END SUBROUTINE fill_tg

!-----------------------------------------------------------------
      SUBROUTINE reset_fnxt()
      INTEGER :: ic
      DO ic = 1,3
         fnxt(ic)%ccomp = field(ic)%ccomp
      END DO
      CALL GState_update_to(fnxt)
      END SUBROUTINE reset_fnxt

!-----------------------------------------------------------------
      SUBROUTINE check_c(name,a,b,tol)
!
! Relative max-norm difference between two spectral arrays
      CHARACTER(len=*), INTENT(IN) :: name
      COMPLEX(KIND=GP), INTENT(IN) :: a(nz,ny,ista:iend),b(nz,ny,ista:iend)
      REAL(KIND=GP), INTENT(IN) :: tol
      REAL(KIND=GP) :: d(2),dl(2)
      INTEGER :: loc(3)
      dl(1) = MAXVAL(ABS(a-b))
      dl(2) = MAXVAL(ABS(b))
      d = dl
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,d,2,GC_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) CALL verdict(name,d(1)/max(d(2),tiny),tol)
      IF (d(1)/max(d(2),tiny).gt.tol .and. dl(1).eq.d(1)) THEN
         loc = MAXLOC(ABS(a-b))
         loc(3) = loc(3)+ista-1
         PRINT '(A,I3,A,3I5,A,2ES12.4,A,2ES12.4)','    rank',myrank,   &
            ' max diff at (k,j,i)=',loc,'  dev=',a(loc(1),loc(2),loc(3)), &
            '  ref=',b(loc(1),loc(2),loc(3))
      ENDIF
      END SUBROUTINE check_c

      SUBROUTINE check_r(name,a,b,tol)
      CHARACTER(len=*), INTENT(IN) :: name
      REAL(KIND=GP), INTENT(IN) :: a(nx,ny,ksta:kend),b(nx,ny,ksta:kend)
      REAL(KIND=GP), INTENT(IN) :: tol
      REAL(KIND=GP) :: d(2)
      d(1) = MAXVAL(ABS(a-b))
      d(2) = MAXVAL(ABS(b))
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,d,2,GC_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) CALL verdict(name,d(1)/max(d(2),tiny),tol)
      END SUBROUTINE check_r

      SUBROUTINE verdict(name,rel,tol)
      CHARACTER(len=*), INTENT(IN) :: name
      REAL(KIND=GP), INTENT(IN) :: rel,tol
      IF (rel.le.tol) THEN
         PRINT '(A,T44,A,ES10.2)',name,'PASS  rel ',rel
      ELSE
         PRINT '(A,T44,A,ES10.2)',name,'FAIL  rel ',rel
         nfail = nfail+1
      ENDIF
      END SUBROUTINE verdict

      END PROGRAM gpu_miniapp
