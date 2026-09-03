!=================================================================
! MODULE gdevice
!
! Device management for OpenMP offload builds (GHOST_GPU): device
! selection per MPI task, synchronization with the vendor runtime,
! creation of device copies of host arrays, transfers between the
! two copies, and the run time switch gdev_active that tells the
! kernels and the FFT layer whether field arrays are to be worked
! on the device or on the host.
!
! Device copies are not created with "target enter data": for a
! derived type component as list item flang maps the whole chain
! of enclosing descriptors (including the CLASS descriptor of a
! "this" dummy, which lives on the caller's stack), and those
! entries later collide with the descriptors kernels map. Instead
! the device memory is allocated with the vendor allocator and
! registered in the OpenMP device table with
! omp_target_associate_ptr on the raw data address, so only the
! data ever sits in the table. Memory from omp_target_alloc is
! not used because GPU-aware MPI moves it about 30 times slower.
!
! Everything here works on addresses and byte counts and has no
! dependency on the rest of GHOST; the typed wrappers are in the
! module gmem. In host builds every routine is a no-op and
! gdev_active is always false.
!
! 2026 Pablo D. Mininni
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!=================================================================

#if defined(GPU_NVIDIA)
#  define C_DEV_SETDEVICE "cudaSetDevice"
#  define C_DEV_SYNC      "cudaDeviceSynchronize"
#  define C_DEV_MALLOC    "cudaMalloc"
#  define C_DEV_FREE      "cudaFree"
#elif defined(GPU_AMD)
#  define C_DEV_SETDEVICE "hipSetDevice"
#  define C_DEV_SYNC      "hipDeviceSynchronize"
#  define C_DEV_MALLOC    "hipMalloc"
#  define C_DEV_FREE      "hipFree"
#elif defined(GHOST_GPU)
#  error 'MODULE GDEVICE: GHOST_GPU needs GPU_NVIDIA or GPU_AMD'
#endif

      MODULE gdevice
      USE ISO_C_BINDING
!$    USE omp_lib
      IMPLICIT NONE

      INTEGER, SAVE :: ndev        = 0       ! devices in this node
      INTEGER, SAVE :: targetdev   = -1      ! device of this task
      INTEGER, SAVE :: hostdev     = -1      ! the host as a device
      LOGICAL, SAVE :: gdev_active = .FALSE. ! fields are worked on the device

! Registry of the device copies: host address, device address and
! size of every array created with gdev_alloc. The registry is what
! the transfers use to find the device copy (omp_get_mapped_ptr does
! not return the pointers registered with omp_target_associate_ptr
! in all runtimes). A few dozen arrays at most, so a linear search.
      INTEGER, PARAMETER :: GDEV_MAXREG = 4096
      INTEGER, SAVE      :: nreg = 0
      TYPE(C_PTR)      , SAVE :: reg_hp(GDEV_MAXREG), reg_dp(GDEV_MAXREG)
      INTEGER(C_SIZE_T), SAVE :: reg_nb(GDEV_MAXREG)

#if defined(GHOST_GPU)
      INTERFACE
        INTEGER(C_INT) FUNCTION dev_setdevice_c(dev) BIND(C,NAME=C_DEV_SETDEVICE)
          IMPORT :: C_INT
          INTEGER(C_INT), VALUE :: dev
        END FUNCTION dev_setdevice_c
        INTEGER(C_INT) FUNCTION dev_sync_c() BIND(C,NAME=C_DEV_SYNC)
          IMPORT :: C_INT
        END FUNCTION dev_sync_c
        INTEGER(C_INT) FUNCTION dev_malloc_c(dp,sz) BIND(C,NAME=C_DEV_MALLOC)
          IMPORT :: C_INT,C_PTR,C_SIZE_T
          TYPE(C_PTR)      , INTENT(OUT) :: dp
          INTEGER(C_SIZE_T), VALUE       :: sz
        END FUNCTION dev_malloc_c
        INTEGER(C_INT) FUNCTION dev_free_c(dp) BIND(C,NAME=C_DEV_FREE)
          IMPORT :: C_INT,C_PTR
          TYPE(C_PTR), VALUE :: dp
        END FUNCTION dev_free_c
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
        INTEGER(C_INT) FUNCTION omp_tgt_memcpy(dst,src,length,dstoff,srcoff,dstdev,srcdev) &
                                BIND(C,NAME="omp_target_memcpy")
          IMPORT :: C_INT,C_PTR,C_SIZE_T
          TYPE(C_PTR)      , VALUE :: dst,src
          INTEGER(C_SIZE_T), VALUE :: length,dstoff,srcoff
          INTEGER(C_INT)   , VALUE :: dstdev,srcdev
        END FUNCTION omp_tgt_memcpy
      END INTERFACE
#endif

      CONTAINS

!*****************************************************************
      SUBROUTINE device_init(myrank)
!-----------------------------------------------------------------
!
! Binds this MPI task to one device, targetdev = myrank modulo the
! number of devices in the node, for OpenMP and for the vendor
! runtime (used by the FFT library and by the allocator). GHOST
! assumes the number of MPI tasks per node is a multiple of the
! number of devices per node; the user must ensure this.
!
! Parameters
!     myrank : MPI rank of this task [IN]
!-----------------------------------------------------------------
      INTEGER, INTENT(IN) :: myrank

#if defined(GHOST_GPU)
!$    ndev    = omp_get_num_devices()
!$    hostdev = omp_get_initial_device()
      IF (ndev.le.0) THEN
         PRINT *,'device_init: no devices found in the node'
         STOP
      ENDIF
      targetdev = MODULO(myrank,ndev)
!$    CALL omp_set_default_device(targetdev)
      IF (dev_setdevice_c(INT(targetdev,C_INT)).ne.0) THEN
         PRINT *,'device_init: the vendor runtime could not select device ',targetdev
         STOP
      ENDIF
#endif
      END SUBROUTINE device_init

!*****************************************************************
      SUBROUTINE device_sync()
!-----------------------------------------------------------------
!
! Waits for all work on the device, including work started by the
! vendor FFT library outside the OpenMP runtime.
!-----------------------------------------------------------------
#if defined(GHOST_GPU)
      INTEGER(C_INT) :: ier
      ier = dev_sync_c()
      IF (ier.ne.0) PRINT *,'device_sync: vendor runtime error ',ier
#endif
      END SUBROUTINE device_sync

!*****************************************************************
      SUBROUTINE gdev_alloc(hp,nbytes)
!-----------------------------------------------------------------
!
! Creates the device copy of a host array and associates it with
! the host address, so that kernels find the array present.
!
! Parameters
!     hp     : host address of the array [IN]
!     nbytes : size in bytes [IN]
!-----------------------------------------------------------------
      TYPE(C_PTR)      , INTENT(IN) :: hp
      INTEGER(C_SIZE_T), INTENT(IN) :: nbytes
#if defined(GHOST_GPU)
      TYPE(C_PTR)    :: dp
      INTEGER(C_INT) :: ier
      IF (nbytes.le.0) RETURN
      ier = dev_malloc_c(dp,nbytes)
      IF (ier.ne.0) THEN
         PRINT *,'gdev_alloc: device allocation of ',nbytes,' bytes failed, error ',ier
         STOP
      ENDIF
      ier = omp_tgt_assoc(hp,dp,nbytes,0_C_SIZE_T,INT(targetdev,C_INT))
      IF (ier.ne.0) THEN
         PRINT *,'gdev_alloc: omp_target_associate_ptr failed, error ',ier
         STOP
      ENDIF
      IF (nreg.ge.GDEV_MAXREG) THEN
         PRINT *,'gdev_alloc: too many device arrays, increase GDEV_MAXREG'
         STOP
      ENDIF
      nreg = nreg+1
      reg_hp(nreg) = hp
      reg_dp(nreg) = dp
      reg_nb(nreg) = nbytes
#endif
      END SUBROUTINE gdev_alloc

!*****************************************************************
      SUBROUTINE gdev_free(hp)
!-----------------------------------------------------------------
!
! Releases the device copy associated with a host address.
!
! Parameters
!     hp : host address of the array [IN]
!-----------------------------------------------------------------
      TYPE(C_PTR), INTENT(IN) :: hp
#if defined(GHOST_GPU)
      INTEGER(C_INT) :: ier
      INTEGER        :: i
      i = gdev_lookup(hp)
      IF (i.eq.0) RETURN
      ier = omp_tgt_disassoc(hp,INT(targetdev,C_INT))
      ier = dev_free_c(reg_dp(i))
      reg_hp(i) = reg_hp(nreg)
      reg_dp(i) = reg_dp(nreg)
      reg_nb(i) = reg_nb(nreg)
      nreg = nreg-1
#endif
      END SUBROUTINE gdev_free

!*****************************************************************
      INTEGER FUNCTION gdev_lookup(hp)
!-----------------------------------------------------------------
!
! Index in the registry of the array at host address hp, 0 if the
! address has no device copy.
!-----------------------------------------------------------------
      TYPE(C_PTR), INTENT(IN) :: hp
      INTEGER :: i
      gdev_lookup = 0
      DO i = 1,nreg
         IF (C_ASSOCIATED(hp,reg_hp(i))) THEN
            gdev_lookup = i
            RETURN
         ENDIF
      END DO
      END FUNCTION gdev_lookup

!*****************************************************************
      SUBROUTINE gdev_update_to(hp,nbytes)
!-----------------------------------------------------------------
!
! Copies the host array to its device copy.
!-----------------------------------------------------------------
      TYPE(C_PTR)      , INTENT(IN) :: hp
      INTEGER(C_SIZE_T), INTENT(IN) :: nbytes
#if defined(GHOST_GPU)
      INTEGER(C_INT) :: ier
      INTEGER        :: i
      i = gdev_lookup(hp)
      IF (i.eq.0) THEN
         PRINT *,'gdev_update_to: array has no device copy'
         STOP
      ENDIF
      ier = omp_tgt_memcpy(reg_dp(i),hp,nbytes,0_C_SIZE_T,0_C_SIZE_T, &
                           INT(targetdev,C_INT),INT(hostdev,C_INT))
      IF (ier.ne.0) PRINT *,'gdev_update_to: omp_target_memcpy error ',ier
#endif
      END SUBROUTINE gdev_update_to

!*****************************************************************
      SUBROUTINE gdev_update_from(hp,nbytes)
!-----------------------------------------------------------------
!
! Copies the device copy of an array back to the host array.
!-----------------------------------------------------------------
      TYPE(C_PTR)      , INTENT(IN) :: hp
      INTEGER(C_SIZE_T), INTENT(IN) :: nbytes
#if defined(GHOST_GPU)
      INTEGER(C_INT) :: ier
      INTEGER        :: i
      i = gdev_lookup(hp)
      IF (i.eq.0) THEN
         PRINT *,'gdev_update_from: array has no device copy'
         STOP
      ENDIF
      ier = omp_tgt_memcpy(hp,reg_dp(i),nbytes,0_C_SIZE_T,0_C_SIZE_T, &
                           INT(hostdev,C_INT),INT(targetdev,C_INT))
      IF (ier.ne.0) PRINT *,'gdev_update_from: omp_target_memcpy error ',ier
#endif
      END SUBROUTINE gdev_update_from

!*****************************************************************
      SUBROUTINE gdev_copy(hdst,hsrc,nbytes)
!-----------------------------------------------------------------
!
! Copies between the device copies of two host arrays.
!-----------------------------------------------------------------
      TYPE(C_PTR)      , INTENT(IN) :: hdst,hsrc
      INTEGER(C_SIZE_T), INTENT(IN) :: nbytes
#if defined(GHOST_GPU)
      INTEGER(C_INT) :: ier
      INTEGER        :: id,is
      id = gdev_lookup(hdst)
      is = gdev_lookup(hsrc)
      IF ((id.eq.0).or.(is.eq.0)) THEN
         PRINT *,'gdev_copy: array has no device copy'
         STOP
      ENDIF
      ier = omp_tgt_memcpy(reg_dp(id),reg_dp(is),nbytes,0_C_SIZE_T,   &
                           0_C_SIZE_T,INT(targetdev,C_INT),INT(targetdev,C_INT))
      IF (ier.ne.0) PRINT *,'gdev_copy: omp_target_memcpy error ',ier
#endif
      END SUBROUTINE gdev_copy

!*****************************************************************
      LOGICAL FUNCTION gdev_present(hp)
!-----------------------------------------------------------------
!
! True if the host address has a device copy (always true in host
! builds, where host and "device" copy are the same array).
!-----------------------------------------------------------------
      TYPE(C_PTR), INTENT(IN) :: hp
#if defined(GHOST_GPU)
      gdev_present = (gdev_lookup(hp).ne.0)
#else
      gdev_present = .TRUE.
#endif
      END FUNCTION gdev_present

      END MODULE gdevice
