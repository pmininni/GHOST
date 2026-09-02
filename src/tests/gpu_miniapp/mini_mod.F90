!=================================================================
! mini_mod.F90
! Basic modules of the GHOST GPU mini-app: precision, MPI ranges,
! grid, wavenumbers, field states, and the workspace pool. The
! types and allocation choke points mirror GHOST (gstate_mod.f90,
! workspace.f90, pseudospec3D_mod.f90) so that what is learned
! here about device residency transfers directly.
!=================================================================
#include "gpu_defs.h"

!=================================================================
  MODULE mini_prec
      USE ISO_C_BINDING
#if defined(GDOUBLE_PRECISION)
      INTEGER, PARAMETER :: GP = KIND(0.0D0)
#else
      INTEGER, PARAMETER :: GP = KIND(0.0)
#endif
  END MODULE mini_prec

!=================================================================
  MODULE mini_mpivars
      USE mpi_f08
      INTEGER, SAVE :: ista,iend
      INTEGER, SAVE :: ksta,kend
      INTEGER, SAVE :: nprocs,myrank
      INTEGER, SAVE :: provided
      INTEGER, SAVE :: ierr
      INTEGER, SAVE :: ndev = 0, targetdev = -1
      TYPE(MPI_Datatype), SAVE :: GC_REAL, GC_COMPLEX
    CONTAINS
      SUBROUTINE range(n1,n2,nprocs,irank,ista,iend)
      INTEGER, INTENT(IN)  :: n1,n2,nprocs,irank
      INTEGER, INTENT(OUT) :: ista,iend
      INTEGER :: iwork1,iwork2
      iwork1 = (n2-n1+1)/nprocs
      iwork2 = MOD(n2-n1+1,nprocs)
      ista = irank*iwork1+n1+MIN(irank,iwork2)
      iend = ista+iwork1-1
      IF (iwork2.gt.irank) iend = iend+1
      END SUBROUTINE range
  END MODULE mini_mpivars

!=================================================================
  MODULE mini_grid
      INTEGER, SAVE :: nx,ny,nz
  END MODULE mini_grid

!=================================================================
  MODULE mini_device
!
! Device management: device selection, synchronization with the
! vendor runtime (needed after FFT library calls), a query for the
! presence of a host address in the OpenMP device table, and the
! helpers that create device copies of host arrays.
      USE ISO_C_BINDING
      USE mini_prec
      USE mini_mpivars
      IMPLICIT NONE
#if defined(GHOST_GPU)
      INTERFACE
        INTEGER(C_INT) FUNCTION dev_sync_c() BIND(C,NAME=C_DEV_SYNC)
          IMPORT :: C_INT
        END FUNCTION dev_sync_c
        INTEGER(C_INT) FUNCTION dev_setdevice_c(dev) BIND(C,NAME=C_DEV_SETDEVICE)
          IMPORT :: C_INT
          INTEGER(C_INT), VALUE :: dev
        END FUNCTION dev_setdevice_c
        INTEGER(C_INT) FUNCTION omp_tgt_is_present(ptr,dev) &
                                BIND(C,NAME="omp_target_is_present")
          IMPORT :: C_INT,C_PTR
          TYPE(C_PTR)   , VALUE :: ptr
          INTEGER(C_INT), VALUE :: dev
        END FUNCTION omp_tgt_is_present
        TYPE(C_PTR) FUNCTION omp_tgt_alloc(sz,dev) BIND(C,NAME="omp_target_alloc")
          IMPORT :: C_INT,C_PTR,C_SIZE_T
          INTEGER(C_SIZE_T), VALUE :: sz
          INTEGER(C_INT)   , VALUE :: dev
        END FUNCTION omp_tgt_alloc
        SUBROUTINE omp_tgt_free(p,dev) BIND(C,NAME="omp_target_free")
          IMPORT :: C_INT,C_PTR
          TYPE(C_PTR)   , VALUE :: p
          INTEGER(C_INT), VALUE :: dev
        END SUBROUTINE omp_tgt_free
        INTEGER(C_INT) FUNCTION omp_tgt_assoc(hp,dp,sz,off,dev) &
                                BIND(C,NAME="omp_target_associate_ptr")
          IMPORT :: C_INT,C_PTR,C_SIZE_T
          TYPE(C_PTR)      , VALUE :: hp,dp
          INTEGER(C_SIZE_T), VALUE :: sz,off
          INTEGER(C_INT)   , VALUE :: dev
        END FUNCTION omp_tgt_assoc
        INTEGER(C_INT) FUNCTION omp_tgt_disassoc(hp,dev) &
                                BIND(C,NAME="omp_target_disassociate_ptr")
          IMPORT :: C_INT,C_PTR
          TYPE(C_PTR)   , VALUE :: hp
          INTEGER(C_INT), VALUE :: dev
        END FUNCTION omp_tgt_disassoc
        TYPE(C_PTR) FUNCTION omp_get_mapped(hp,dev) BIND(C,NAME="omp_get_mapped_ptr")
          IMPORT :: C_INT,C_PTR
          TYPE(C_PTR)   , VALUE :: hp
          INTEGER(C_INT), VALUE :: dev
        END FUNCTION omp_get_mapped
#if defined(USE_VENDOR_ALLOC)
! Vendor allocator (hipMalloc / cudaMalloc): same C signature
        INTEGER(C_INT) FUNCTION vendor_malloc(dp,sz) BIND(C,NAME=C_DEV_MALLOC)
          IMPORT :: C_INT,C_PTR,C_SIZE_T
          TYPE(C_PTR)      , INTENT(OUT) :: dp
          INTEGER(C_SIZE_T), VALUE :: sz
        END FUNCTION vendor_malloc
        INTEGER(C_INT) FUNCTION vendor_free(dp) BIND(C,NAME=C_DEV_FREE)
          IMPORT :: C_INT,C_PTR
          TYPE(C_PTR), VALUE :: dp
        END FUNCTION vendor_free
#endif
      END INTERFACE
#endif
    CONTAINS
      SUBROUTINE device_init(myrank)
!$    USE omp_lib
      INTEGER, INTENT(IN) :: myrank
      ndev = 0
      targetdev = -1
#if defined(GHOST_GPU)
!$    ndev = omp_get_num_devices()
      IF (ndev.eq.0) THEN
         PRINT *,'device_init: no devices found'
         STOP
      ENDIF
      targetdev = MODULO(myrank,ndev)
!$    CALL omp_set_default_device(targetdev)
! The vendor runtime (used by the FFT library and by the vendor
! allocator) keeps its own notion of current device; it must be the
! same device as the OpenMP default device.
      IF (dev_setdevice_c(INT(targetdev,C_INT)).ne.0) &
         STOP 'device_init: vendor set-device call failed'
#endif
      END SUBROUTINE device_init

      SUBROUTINE device_sync()
#if defined(GHOST_GPU)
      INTEGER(C_INT) :: ier
      ier = dev_sync_c()
      IF (ier.ne.0) PRINT *,'device_sync: runtime error ',ier
#endif
      END SUBROUTINE device_sync

      LOGICAL FUNCTION device_present(p)
      TYPE(C_PTR), INTENT(IN) :: p
#if defined(GHOST_GPU)
      device_present = (omp_tgt_is_present(p,INT(targetdev,C_INT)).ne.0)
#else
      device_present = .TRUE.
#endif
      END FUNCTION device_present

! ----------------------------------------------------------------
! Device residency helpers. A device copy of a host array is created
! with omp_target_alloc and registered in the device table with
! omp_target_associate_ptr, on the raw data address and byte count.
! This is used instead of "target enter data" at the allocation
! choke points because, with a derived type component as list item,
! flang maps the whole chain of enclosing descriptors (including the
! CLASS(...) descriptor of a "this" dummy, which lives on the
! caller's stack) and those persistent stack entries later collide
! with the transient descriptors that kernels map. The entries made
! here hold only the data, exactly as the design wants.
! ----------------------------------------------------------------
      SUBROUTINE dev_alloc_bytes(hp,nbytes)
      TYPE(C_PTR), INTENT(IN) :: hp
      INTEGER(C_SIZE_T), INTENT(IN) :: nbytes
#if defined(GHOST_GPU)
      TYPE(C_PTR) :: dp
      INTEGER(C_INT) :: ier
#if defined(USE_VENDOR_ALLOC)
      ier = vendor_malloc(dp,nbytes)
      IF (ier.ne.0) STOP 'dev_alloc_bytes: vendor malloc failed'
#else
      dp = omp_tgt_alloc(nbytes,INT(targetdev,C_INT))
      IF (.not.C_ASSOCIATED(dp)) STOP 'dev_alloc_bytes: omp_target_alloc failed'
#endif
      ier = omp_tgt_assoc(hp,dp,nbytes,0_C_SIZE_T,INT(targetdev,C_INT))
      IF (ier.ne.0) STOP 'dev_alloc_bytes: omp_target_associate_ptr failed'
#endif
      END SUBROUTINE dev_alloc_bytes

      SUBROUTINE dev_free_bytes(hp)
      TYPE(C_PTR), INTENT(IN) :: hp
#if defined(GHOST_GPU)
      TYPE(C_PTR) :: dp
      INTEGER(C_INT) :: ier
      dp  = omp_get_mapped(hp,INT(targetdev,C_INT))
      ier = omp_tgt_disassoc(hp,INT(targetdev,C_INT))
#if defined(USE_VENDOR_ALLOC)
      IF (C_ASSOCIATED(dp)) ier = vendor_free(dp)
#else
      IF (C_ASSOCIATED(dp)) CALL omp_tgt_free(dp,INT(targetdev,C_INT))
#endif
#endif
      END SUBROUTINE dev_free_bytes

      SUBROUTINE dev_alloc_c(a,n)
      INTEGER, INTENT(IN) :: n
      COMPLEX(KIND=GP), INTENT(IN), TARGET :: a(n)
      CALL dev_alloc_bytes(C_LOC(a),INT(n,C_SIZE_T)*INT(2*GP,C_SIZE_T))
      END SUBROUTINE dev_alloc_c

      SUBROUTINE dev_free_c(a,n)
      INTEGER, INTENT(IN) :: n
      COMPLEX(KIND=GP), INTENT(IN), TARGET :: a(n)
      CALL dev_free_bytes(C_LOC(a))
      END SUBROUTINE dev_free_c

      SUBROUTINE dev_alloc_r(a,n)
      INTEGER, INTENT(IN) :: n
      REAL(KIND=GP), INTENT(IN), TARGET :: a(n)
      CALL dev_alloc_bytes(C_LOC(a),INT(n,C_SIZE_T)*INT(GP,C_SIZE_T))
      END SUBROUTINE dev_alloc_r

      SUBROUTINE dev_free_r(a,n)
      INTEGER, INTENT(IN) :: n
      REAL(KIND=GP), INTENT(IN), TARGET :: a(n)
      CALL dev_free_bytes(C_LOC(a))
      END SUBROUTINE dev_free_r

      SUBROUTINE dev_update_to_r(a,n)
!
! Host to device copy of a real array of any rank (explicit-shape
! dummy, so only the data range is looked up in the device table)
      INTEGER, INTENT(IN) :: n
      REAL(KIND=GP), INTENT(IN) :: a(n)
#if defined(GHOST_GPU)
!$omp target update to(a)
#endif
      END SUBROUTINE dev_update_to_r
  END MODULE mini_device

!=================================================================
  MODULE mini_xfer
!
! Explicit host<->device updates of whole field-sized arrays. The
! dummies are explicit-shape, so a state component or a pool array
! can be passed directly; the update resolves the address in the
! device table. In a host build these are no-ops.
      USE mini_prec
      USE mini_grid
      USE mini_mpivars
      IMPLICIT NONE
    CONTAINS
      SUBROUTINE cupdate_to(a)
      COMPLEX(KIND=GP), INTENT(IN) :: a(nz,ny,ista:iend)
#if defined(GHOST_GPU)
!$omp target update to(a)
#endif
      END SUBROUTINE cupdate_to

      SUBROUTINE cupdate_from(a)
      COMPLEX(KIND=GP), INTENT(INOUT) :: a(nz,ny,ista:iend)
#if defined(GHOST_GPU)
!$omp target update from(a)
#endif
      END SUBROUTINE cupdate_from

      SUBROUTINE rupdate_to(a)
      REAL(KIND=GP), INTENT(IN) :: a(nx,ny,ksta:kend)
#if defined(GHOST_GPU)
!$omp target update to(a)
#endif
      END SUBROUTINE rupdate_to

      SUBROUTINE rupdate_from(a)
      REAL(KIND=GP), INTENT(INOUT) :: a(nx,ny,ksta:kend)
#if defined(GHOST_GPU)
!$omp target update from(a)
#endif
      END SUBROUTINE rupdate_from
  END MODULE mini_xfer

!=================================================================
  MODULE mini_kes
!
! Wavenumbers for a 2.pi cubic box, as in box_init. The arrays
! are module variables referenced from inside device kernels; they
! are mapped once here and never updated.
      USE mini_prec
      USE mini_grid
      USE mini_mpivars
      USE mini_device
      IMPLICIT NONE
      REAL(KIND=GP), ALLOCATABLE, SAVE :: kx(:),ky(:),kz(:)
      REAL(KIND=GP), ALLOCATABLE, SAVE :: kk2(:,:,:)
      REAL(KIND=GP), SAVE :: kmax,tiny
    CONTAINS
      SUBROUTINE kes_init()
      INTEGER :: i,j,k
      ALLOCATE( kx(nx), ky(ny), kz(nz) )
      ALLOCATE( kk2(nz,ny,ista:iend) )
      DO i = 1,nx/2
         kx(i) = real(i-1,kind=GP)
         kx(i+nx/2) = real(i-nx/2-1,kind=GP)
      END DO
      DO j = 1,ny/2
         ky(j) = real(j-1,kind=GP)
         ky(j+ny/2) = real(j-ny/2-1,kind=GP)
      END DO
      DO k = 1,nz/2
         kz(k) = real(k-1,kind=GP)
         kz(k+nz/2) = real(k-nz/2-1,kind=GP)
      END DO
      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               kk2(k,j,i) = kx(i)**2+ky(j)**2+kz(k)**2
            END DO
         END DO
      END DO
      kmax = (real(nx,kind=GP)/3.0_GP)**2
      tiny = min(1e-5_GP,.1_GP/(real(nx,kind=GP)**2))
      CALL dev_alloc_r(kx ,nx);           CALL dev_update_to_r(kx ,nx)
      CALL dev_alloc_r(ky ,ny);           CALL dev_update_to_r(ky ,ny)
      CALL dev_alloc_r(kz ,nz);           CALL dev_update_to_r(kz ,nz)
      CALL dev_alloc_r(kk2,SIZE(kk2));    CALL dev_update_to_r(kk2,SIZE(kk2))
      END SUBROUTINE kes_init
  END MODULE mini_kes

!=================================================================
  MODULE mini_state
!
! Field state: an array of derived types, one allocatable complex
! component each (GHOST's GStateComp). Device residency is decided
! at the single allocation point: the component is allocated on
! the host and an uninitialized device copy is created and
! associated with it. Nothing else in the code maps state arrays.
      USE mini_prec
      USE mini_grid
      USE mini_mpivars
      USE mini_device
      USE mini_xfer
      IMPLICIT NONE
      TYPE, PUBLIC :: GStateComp
        COMPLEX(KIND=GP), ALLOCATABLE :: ccomp(:,:,:)
      END TYPE GStateComp
    CONTAINS
      SUBROUTINE GState_alloc(state,nc)
      TYPE(GStateComp), ALLOCATABLE, INTENT(INOUT) :: state(:)
      INTEGER, INTENT(IN) :: nc
      INTEGER :: i
      ALLOCATE( state(nc) )
      DO i = 1,nc
         ALLOCATE( state(i)%ccomp(nz,ny,ista:iend) )
         state(i)%ccomp = 0.0_GP
         CALL dev_alloc_c(state(i)%ccomp,SIZE(state(i)%ccomp))
      END DO
      END SUBROUTINE GState_alloc

      SUBROUTINE GState_dealloc(state)
      TYPE(GStateComp), ALLOCATABLE, INTENT(INOUT) :: state(:)
      INTEGER :: i
      IF (.not.ALLOCATED(state)) RETURN
      DO i = 1,SIZE(state)
         IF (ALLOCATED(state(i)%ccomp)) THEN
            CALL dev_free_c(state(i)%ccomp,SIZE(state(i)%ccomp))
            DEALLOCATE( state(i)%ccomp )
         ENDIF
      END DO
      DEALLOCATE( state )
      END SUBROUTINE GState_dealloc

      SUBROUTINE GState_update_to(state)
      TYPE(GStateComp), INTENT(IN) :: state(:)
      INTEGER :: i
      DO i = 1,SIZE(state)
         CALL cupdate_to(state(i)%ccomp)
      END DO
      END SUBROUTINE GState_update_to

      SUBROUTINE GState_update_from(state)
      TYPE(GStateComp), INTENT(INOUT) :: state(:)
      INTEGER :: i
      DO i = 1,SIZE(state)
         CALL cupdate_from(state(i)%ccomp)
      END DO
      END SUBROUTINE GState_update_from
  END MODULE mini_state

!=================================================================
  MODULE mini_workspace
!
! Workspace pool (GHOST's GWorkspace): preallocated real and
! complex arrays handed out as pointers. As with the states, the
! device copies are created once at initialize_pool and kernels
! find them through the device table when they receive a pointer.
      USE mini_prec
      USE mini_grid
      USE mini_mpivars
      USE mini_device
      IMPLICIT NONE
      PRIVATE

      TYPE :: RealEntry
        LOGICAL :: is_free = .TRUE.
        REAL(KIND=GP), ALLOCATABLE :: array(:,:,:)
      END TYPE RealEntry
      TYPE :: ComplexEntry
        LOGICAL :: is_free = .TRUE.
        COMPLEX(KIND=GP), ALLOCATABLE :: array(:,:,:)
      END TYPE ComplexEntry

      TYPE, PUBLIC :: GWorkspace
        TYPE(RealEntry)   , ALLOCATABLE :: real_entries_(:)
        TYPE(ComplexEntry), ALLOCATABLE :: complex_entries_(:)
      CONTAINS
        PROCEDURE, PUBLIC :: initialize_pool
        PROCEDURE, PUBLIC :: get_real_tmp    , get_complex_tmp
        PROCEDURE, PUBLIC :: free_real_tmp   , free_complex_tmp
        PROCEDURE, PUBLIC :: destroy_pool
      END TYPE GWorkspace

    CONTAINS
      SUBROUTINE initialize_pool(this,num_real,num_complex)
      CLASS(GWorkspace), INTENT(INOUT), TARGET :: this
      INTEGER, INTENT(IN) :: num_real,num_complex
      INTEGER :: i
      ALLOCATE( this%real_entries_(num_real) )
      ALLOCATE( this%complex_entries_(num_complex) )
      DO i = 1,num_real
         ALLOCATE( this%real_entries_(i)%array(nx,ny,ksta:kend) )
         this%real_entries_(i)%array = 0.0_GP
         this%real_entries_(i)%is_free = .TRUE.
         CALL dev_alloc_r(this%real_entries_(i)%array,SIZE(this%real_entries_(i)%array))
      END DO
      DO i = 1,num_complex
         ALLOCATE( this%complex_entries_(i)%array(nz,ny,ista:iend) )
         this%complex_entries_(i)%array = 0.0_GP
         this%complex_entries_(i)%is_free = .TRUE.
         CALL dev_alloc_c(this%complex_entries_(i)%array,SIZE(this%complex_entries_(i)%array))
      END DO
      END SUBROUTINE initialize_pool

      SUBROUTINE destroy_pool(this)
      CLASS(GWorkspace), INTENT(INOUT), TARGET :: this
      INTEGER :: i
      DO i = 1,SIZE(this%real_entries_)
         CALL dev_free_r(this%real_entries_(i)%array,SIZE(this%real_entries_(i)%array))
         DEALLOCATE( this%real_entries_(i)%array )
      END DO
      DO i = 1,SIZE(this%complex_entries_)
         CALL dev_free_c(this%complex_entries_(i)%array,SIZE(this%complex_entries_(i)%array))
         DEALLOCATE( this%complex_entries_(i)%array )
      END DO
      DEALLOCATE( this%real_entries_, this%complex_entries_ )
      END SUBROUTINE destroy_pool

      SUBROUTINE get_real_tmp(this,ptr,bret)
      CLASS(GWorkspace), INTENT(INOUT), TARGET :: this
      REAL(KIND=GP), POINTER, INTENT(OUT) :: ptr(:,:,:)
      LOGICAL, INTENT(OUT) :: bret
      INTEGER :: i
      bret = .FALSE.
      ptr => NULL()
      DO i = 1,SIZE(this%real_entries_)
         IF (this%real_entries_(i)%is_free) THEN
            this%real_entries_(i)%is_free = .FALSE.
            ptr => this%real_entries_(i)%array
            bret = .TRUE.
            RETURN
         ENDIF
      END DO
      STOP 'get_real_tmp: pool exhausted'
      END SUBROUTINE get_real_tmp

      SUBROUTINE get_complex_tmp(this,ptr,bret)
      CLASS(GWorkspace), INTENT(INOUT), TARGET :: this
      COMPLEX(KIND=GP), POINTER, INTENT(OUT) :: ptr(:,:,:)
      LOGICAL, INTENT(OUT) :: bret
      INTEGER :: i
      bret = .FALSE.
      ptr => NULL()
      DO i = 1,SIZE(this%complex_entries_)
         IF (this%complex_entries_(i)%is_free) THEN
            this%complex_entries_(i)%is_free = .FALSE.
            ptr => this%complex_entries_(i)%array
            bret = .TRUE.
            RETURN
         ENDIF
      END DO
      STOP 'get_complex_tmp: pool exhausted'
      END SUBROUTINE get_complex_tmp

      SUBROUTINE free_real_tmp(this,ptr)
      CLASS(GWorkspace), INTENT(INOUT), TARGET :: this
      REAL(KIND=GP), POINTER, INTENT(INOUT) :: ptr(:,:,:)
      INTEGER :: i
      DO i = 1,SIZE(this%real_entries_)
         IF (ASSOCIATED(ptr,this%real_entries_(i)%array)) THEN
            this%real_entries_(i)%is_free = .TRUE.
            ptr => NULL()
            RETURN
         ENDIF
      END DO
      STOP 'free_real_tmp: pointer not from pool'
      END SUBROUTINE free_real_tmp

      SUBROUTINE free_complex_tmp(this,ptr)
      CLASS(GWorkspace), INTENT(INOUT), TARGET :: this
      COMPLEX(KIND=GP), POINTER, INTENT(INOUT) :: ptr(:,:,:)
      INTEGER :: i
      DO i = 1,SIZE(this%complex_entries_)
         IF (ASSOCIATED(ptr,this%complex_entries_(i)%array)) THEN
            this%complex_entries_(i)%is_free = .TRUE.
            ptr => NULL()
            RETURN
         ENDIF
      END DO
      STOP 'free_complex_tmp: pointer not from pool'
      END SUBROUTINE free_complex_tmp
  END MODULE mini_workspace
