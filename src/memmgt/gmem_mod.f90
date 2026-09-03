!=================================================================
! MODULE gmem
!
! Typed allocator for the field-sized arrays of GHOST. Every array
! that a kernel may touch (states, workspace pool, wavenumbers,
! FFT buffers) is allocated through galloc and released through
! gfree, so that the decision of where an array lives is taken in
! one place: on the host in host builds, on the host with a device
! copy in offload builds. The remaining routines move data between
! the two copies:
!
!   galloc(a,...)     : ALLOCATE plus device copy (rank 1 or 3,
!                       real or complex, allocatable or pointer)
!   galloc_host(a,...): ALLOCATE without device copy (host-only arrays)
!   gfree(a)          : release both copies
!   gupdate_to(a)     : host copy -> device copy
!   gupdate_from(a)   : device copy -> host copy
!   gcopy(dst,src)    : copy between arrays; in offload builds the
!                       device copies are copied and the host copies
!                       are left untouched, in host builds it is an
!                       assignment
!
! In host builds galloc/gfree reduce to ALLOCATE/DEALLOCATE and the
! updates are no-ops. Rank 3 arrays take the bounds (n1,n2,l3:u3),
! matching the (nz,ny,ista:iend) and (nx,ny,ksta:kend) layouts.
!
! 2026 Pablo D. Mininni
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!=================================================================

      MODULE gmem
      USE ISO_C_BINDING
      USE fprecision
      USE gdevice
      IMPLICIT NONE

      INTERFACE galloc
        MODULE PROCEDURE galloc_r1, galloc_r3, galloc_c1, galloc_c3, &
                         galloc_i1, galloc_i2, galloc_r2
      END INTERFACE galloc
      INTERFACE galloc_ptr
        MODULE PROCEDURE galloc_r3p
      END INTERFACE galloc_ptr
      INTERFACE galloc_host
        MODULE PROCEDURE galloc_host_r3, galloc_host_c3
      END INTERFACE galloc_host
      INTERFACE gfree
        MODULE PROCEDURE gfree_r1, gfree_r3, gfree_c1, gfree_c3, &
                         gfree_i1, gfree_i2, gfree_r2
      END INTERFACE gfree
      INTERFACE gfree_ptr
        MODULE PROCEDURE gfree_r3p
      END INTERFACE gfree_ptr
      INTERFACE gupdate_to
        MODULE PROCEDURE gupdate_to_r1, gupdate_to_r3, gupdate_to_c1, gupdate_to_c3, &
                         gupdate_to_i1, gupdate_to_i2, gupdate_to_r2
      END INTERFACE gupdate_to
      INTERFACE gupdate_from
        MODULE PROCEDURE gupdate_from_r1, gupdate_from_r3, gupdate_from_c1, gupdate_from_c3, &
                         gupdate_from_i1, gupdate_from_i2, gupdate_from_r2
      END INTERFACE gupdate_from
      INTERFACE gcopy
        MODULE PROCEDURE gcopy_r3, gcopy_c3, gcopy_r1, gcopy_i1, gcopy_r2, gcopy_i2
      END INTERFACE gcopy
      INTERFACE gresize
        MODULE PROCEDURE gresize_r1, gresize_i1, gresize_r2, gresize_i2
      END INTERFACE gresize

      CONTAINS

! ---------------------------------------------------------------
! Allocation
! ---------------------------------------------------------------
      SUBROUTINE galloc_r1(a,n)
      REAL(KIND=GP), ALLOCATABLE, TARGET, INTENT(OUT) :: a(:)
      INTEGER, INTENT(IN) :: n
      ALLOCATE( a(n) )
      CALL gdev_alloc(C_LOC(a),nbytes_r(SIZE(a)))
      END SUBROUTINE galloc_r1

      SUBROUTINE galloc_r3(a,n1,n2,l3,u3)
      REAL(KIND=GP), ALLOCATABLE, TARGET, INTENT(OUT) :: a(:,:,:)
      INTEGER, INTENT(IN) :: n1,n2,l3,u3
      ALLOCATE( a(n1,n2,l3:u3) )
      CALL gdev_alloc(C_LOC(a),nbytes_r(SIZE(a)))
      END SUBROUTINE galloc_r3

      SUBROUTINE galloc_r3p(a,n1,n2,l3,u3)
      REAL(KIND=GP), POINTER, INTENT(OUT) :: a(:,:,:)
      INTEGER, INTENT(IN) :: n1,n2,l3,u3
      ALLOCATE( a(n1,n2,l3:u3) )
      CALL gdev_alloc(C_LOC(a),nbytes_r(SIZE(a)))
      END SUBROUTINE galloc_r3p

      SUBROUTINE galloc_c1(a,n)
      COMPLEX(KIND=GP), ALLOCATABLE, TARGET, INTENT(OUT) :: a(:)
      INTEGER, INTENT(IN) :: n
      ALLOCATE( a(n) )
      CALL gdev_alloc(C_LOC(a),nbytes_c(SIZE(a)))
      END SUBROUTINE galloc_c1

      SUBROUTINE galloc_c3(a,n1,n2,l3,u3)
      COMPLEX(KIND=GP), ALLOCATABLE, TARGET, INTENT(OUT) :: a(:,:,:)
      INTEGER, INTENT(IN) :: n1,n2,l3,u3
      ALLOCATE( a(n1,n2,l3:u3) )
      CALL gdev_alloc(C_LOC(a),nbytes_c(SIZE(a)))
      END SUBROUTINE galloc_c3

! ---------------------------------------------------------------
! Host-only allocation: arrays that are only ever used on the host
! (e.g. the temporaries of the diagnostics, which run on the host
! copies of the fields). No device copy is created; gfree releases
! them like any other array.
! ---------------------------------------------------------------
      SUBROUTINE galloc_host_r3(a,n1,n2,l3,u3)
      REAL(KIND=GP), ALLOCATABLE, TARGET, INTENT(OUT) :: a(:,:,:)
      INTEGER, INTENT(IN) :: n1,n2,l3,u3
      ALLOCATE( a(n1,n2,l3:u3) )
      END SUBROUTINE galloc_host_r3

      SUBROUTINE galloc_host_c3(a,n1,n2,l3,u3)
      COMPLEX(KIND=GP), ALLOCATABLE, TARGET, INTENT(OUT) :: a(:,:,:)
      INTEGER, INTENT(IN) :: n1,n2,l3,u3
      ALLOCATE( a(n1,n2,l3:u3) )
      END SUBROUTINE galloc_host_c3

! ---------------------------------------------------------------
! Release
! ---------------------------------------------------------------
      SUBROUTINE gfree_r1(a)
      REAL(KIND=GP), ALLOCATABLE, TARGET, INTENT(INOUT) :: a(:)
      IF (.not.ALLOCATED(a)) RETURN
      CALL gdev_free(C_LOC(a))
      DEALLOCATE( a )
      END SUBROUTINE gfree_r1

      SUBROUTINE gfree_r3(a)
      REAL(KIND=GP), ALLOCATABLE, TARGET, INTENT(INOUT) :: a(:,:,:)
      IF (.not.ALLOCATED(a)) RETURN
      CALL gdev_free(C_LOC(a))
      DEALLOCATE( a )
      END SUBROUTINE gfree_r3

      SUBROUTINE gfree_r3p(a)
      REAL(KIND=GP), POINTER, INTENT(INOUT) :: a(:,:,:)
      IF (.not.ASSOCIATED(a)) RETURN
      CALL gdev_free(C_LOC(a))
      DEALLOCATE( a )
      END SUBROUTINE gfree_r3p

      SUBROUTINE gfree_c1(a)
      COMPLEX(KIND=GP), ALLOCATABLE, TARGET, INTENT(INOUT) :: a(:)
      IF (.not.ALLOCATED(a)) RETURN
      CALL gdev_free(C_LOC(a))
      DEALLOCATE( a )
      END SUBROUTINE gfree_c1

      SUBROUTINE gfree_c3(a)
      COMPLEX(KIND=GP), ALLOCATABLE, TARGET, INTENT(INOUT) :: a(:,:,:)
      IF (.not.ALLOCATED(a)) RETURN
      CALL gdev_free(C_LOC(a))
      DEALLOCATE( a )
      END SUBROUTINE gfree_c3

! ---------------------------------------------------------------
! Host -> device
! ---------------------------------------------------------------
      SUBROUTINE gupdate_to_r1(a)
      REAL(KIND=GP), INTENT(IN), TARGET, CONTIGUOUS :: a(:)
      CALL gdev_update_to(C_LOC(a),nbytes_r(SIZE(a)))
      END SUBROUTINE gupdate_to_r1

      SUBROUTINE gupdate_to_r3(a)
      REAL(KIND=GP), INTENT(IN), TARGET, CONTIGUOUS :: a(:,:,:)
      CALL gdev_update_to(C_LOC(a),nbytes_r(SIZE(a)))
      END SUBROUTINE gupdate_to_r3

      SUBROUTINE gupdate_to_c1(a)
      COMPLEX(KIND=GP), INTENT(IN), TARGET, CONTIGUOUS :: a(:)
      CALL gdev_update_to(C_LOC(a),nbytes_c(SIZE(a)))
      END SUBROUTINE gupdate_to_c1

      SUBROUTINE gupdate_to_c3(a)
      COMPLEX(KIND=GP), INTENT(IN), TARGET, CONTIGUOUS :: a(:,:,:)
      CALL gdev_update_to(C_LOC(a),nbytes_c(SIZE(a)))
      END SUBROUTINE gupdate_to_c3

! ---------------------------------------------------------------
! Device -> host
! ---------------------------------------------------------------
      SUBROUTINE gupdate_from_r1(a)
      REAL(KIND=GP), INTENT(INOUT), TARGET, CONTIGUOUS :: a(:)
      CALL gdev_update_from(C_LOC(a),nbytes_r(SIZE(a)))
      END SUBROUTINE gupdate_from_r1

      SUBROUTINE gupdate_from_r3(a)
      REAL(KIND=GP), INTENT(INOUT), TARGET, CONTIGUOUS :: a(:,:,:)
      CALL gdev_update_from(C_LOC(a),nbytes_r(SIZE(a)))
      END SUBROUTINE gupdate_from_r3

      SUBROUTINE gupdate_from_c1(a)
      COMPLEX(KIND=GP), INTENT(INOUT), TARGET, CONTIGUOUS :: a(:)
      CALL gdev_update_from(C_LOC(a),nbytes_c(SIZE(a)))
      END SUBROUTINE gupdate_from_c1

      SUBROUTINE gupdate_from_c3(a)
      COMPLEX(KIND=GP), INTENT(INOUT), TARGET, CONTIGUOUS :: a(:,:,:)
      CALL gdev_update_from(C_LOC(a),nbytes_c(SIZE(a)))
      END SUBROUTINE gupdate_from_c3

! ---------------------------------------------------------------
! Copies between arrays (device copies in offload builds)
! ---------------------------------------------------------------
      SUBROUTINE gcopy_r3(dst,src)
      REAL(KIND=GP), INTENT(INOUT), TARGET, CONTIGUOUS :: dst(:,:,:)
      REAL(KIND=GP), INTENT(IN)   , TARGET, CONTIGUOUS :: src(:,:,:)
#if defined(GHOST_GPU)
      CALL gdev_copy(C_LOC(dst),C_LOC(src),nbytes_r(SIZE(src)))
#else
      dst = src
#endif
      END SUBROUTINE gcopy_r3

      SUBROUTINE gcopy_c3(dst,src)
      COMPLEX(KIND=GP), INTENT(INOUT), TARGET, CONTIGUOUS :: dst(:,:,:)
      COMPLEX(KIND=GP), INTENT(IN)   , TARGET, CONTIGUOUS :: src(:,:,:)
#if defined(GHOST_GPU)
      CALL gdev_copy(C_LOC(dst),C_LOC(src),nbytes_c(SIZE(src)))
#else
      dst = src
#endif
      END SUBROUTINE gcopy_c3

! ---------------------------------------------------------------
! Sizes in bytes
! ---------------------------------------------------------------
!
! Integer and rank-2 arrays (particle indices, per-particle vectors,
! communication buffers of the particle classes)
!
      SUBROUTINE galloc_i1(a,n)
      INTEGER, ALLOCATABLE, TARGET, INTENT(OUT) :: a(:)
      INTEGER, INTENT(IN) :: n
      ALLOCATE( a(n) )
      CALL gdev_alloc(C_LOC(a),nbytes_i(SIZE(a)))
      END SUBROUTINE galloc_i1

      SUBROUTINE galloc_i2(a,n1,n2)
      INTEGER, ALLOCATABLE, TARGET, INTENT(OUT) :: a(:,:)
      INTEGER, INTENT(IN) :: n1,n2
      ALLOCATE( a(n1,n2) )
      CALL gdev_alloc(C_LOC(a),nbytes_i(SIZE(a)))
      END SUBROUTINE galloc_i2

      SUBROUTINE galloc_r2(a,n1,n2)
      REAL(KIND=GP), ALLOCATABLE, TARGET, INTENT(OUT) :: a(:,:)
      INTEGER, INTENT(IN) :: n1,n2
      ALLOCATE( a(n1,n2) )
      CALL gdev_alloc(C_LOC(a),nbytes_r(SIZE(a)))
      END SUBROUTINE galloc_r2

      SUBROUTINE gfree_i1(a)
      INTEGER, ALLOCATABLE, TARGET, INTENT(INOUT) :: a(:)
      IF (.not.ALLOCATED(a)) RETURN
      CALL gdev_free(C_LOC(a))
      DEALLOCATE( a )
      END SUBROUTINE gfree_i1

      SUBROUTINE gfree_i2(a)
      INTEGER, ALLOCATABLE, TARGET, INTENT(INOUT) :: a(:,:)
      IF (.not.ALLOCATED(a)) RETURN
      CALL gdev_free(C_LOC(a))
      DEALLOCATE( a )
      END SUBROUTINE gfree_i2

      SUBROUTINE gfree_r2(a)
      REAL(KIND=GP), ALLOCATABLE, TARGET, INTENT(INOUT) :: a(:,:)
      IF (.not.ALLOCATED(a)) RETURN
      CALL gdev_free(C_LOC(a))
      DEALLOCATE( a )
      END SUBROUTINE gfree_r2

      SUBROUTINE gupdate_to_i1(a)
      INTEGER, INTENT(IN), TARGET, CONTIGUOUS :: a(:)
      CALL gdev_update_to(C_LOC(a),nbytes_i(SIZE(a)))
      END SUBROUTINE gupdate_to_i1

      SUBROUTINE gupdate_to_i2(a)
      INTEGER, INTENT(IN), TARGET, CONTIGUOUS :: a(:,:)
      CALL gdev_update_to(C_LOC(a),nbytes_i(SIZE(a)))
      END SUBROUTINE gupdate_to_i2

      SUBROUTINE gupdate_to_r2(a)
      REAL(KIND=GP), INTENT(IN), TARGET, CONTIGUOUS :: a(:,:)
      CALL gdev_update_to(C_LOC(a),nbytes_r(SIZE(a)))
      END SUBROUTINE gupdate_to_r2

      SUBROUTINE gupdate_from_i1(a)
      INTEGER, INTENT(INOUT), TARGET, CONTIGUOUS :: a(:)
      CALL gdev_update_from(C_LOC(a),nbytes_i(SIZE(a)))
      END SUBROUTINE gupdate_from_i1

      SUBROUTINE gupdate_from_i2(a)
      INTEGER, INTENT(INOUT), TARGET, CONTIGUOUS :: a(:,:)
      CALL gdev_update_from(C_LOC(a),nbytes_i(SIZE(a)))
      END SUBROUTINE gupdate_from_i2

      SUBROUTINE gupdate_from_r2(a)
      REAL(KIND=GP), INTENT(INOUT), TARGET, CONTIGUOUS :: a(:,:)
      CALL gdev_update_from(C_LOC(a),nbytes_r(SIZE(a)))
      END SUBROUTINE gupdate_from_r2

      SUBROUTINE gcopy_r1(dst,src)
      REAL(KIND=GP), INTENT(INOUT), TARGET, CONTIGUOUS :: dst(:)
      REAL(KIND=GP), INTENT(IN)   , TARGET, CONTIGUOUS :: src(:)
#if defined(GHOST_GPU)
      CALL gdev_copy(C_LOC(dst),C_LOC(src),nbytes_r(SIZE(src)))
#else
      dst = src
#endif
      END SUBROUTINE gcopy_r1

      SUBROUTINE gcopy_i1(dst,src)
      INTEGER, INTENT(INOUT), TARGET, CONTIGUOUS :: dst(:)
      INTEGER, INTENT(IN)   , TARGET, CONTIGUOUS :: src(:)
#if defined(GHOST_GPU)
      CALL gdev_copy(C_LOC(dst),C_LOC(src),nbytes_i(SIZE(src)))
#else
      dst = src
#endif
      END SUBROUTINE gcopy_i1

      SUBROUTINE gcopy_r2(dst,src)
      REAL(KIND=GP), INTENT(INOUT), TARGET, CONTIGUOUS :: dst(:,:)
      REAL(KIND=GP), INTENT(IN)   , TARGET, CONTIGUOUS :: src(:,:)
#if defined(GHOST_GPU)
      CALL gdev_copy(C_LOC(dst),C_LOC(src),nbytes_r(SIZE(src)))
#else
      dst = src
#endif
      END SUBROUTINE gcopy_r2

      SUBROUTINE gcopy_i2(dst,src)
      INTEGER, INTENT(INOUT), TARGET, CONTIGUOUS :: dst(:,:)
      INTEGER, INTENT(IN)   , TARGET, CONTIGUOUS :: src(:,:)
#if defined(GHOST_GPU)
      CALL gdev_copy(C_LOC(dst),C_LOC(src),nbytes_i(SIZE(src)))
#else
      dst = src
#endif
      END SUBROUTINE gcopy_i2

!
! Resizes an array that has a device copy. When keep is true the
! leading MIN(old,new) elements (columns for rank 2) are preserved,
! on the device (the authoritative copy while the time step runs)
! and on the host. Nothing is done if the size does not change.
!
      SUBROUTINE gresize_r1(a,n,keep)
      REAL(KIND=GP), ALLOCATABLE, TARGET, INTENT(INOUT) :: a(:)
      INTEGER, INTENT(IN) :: n
      LOGICAL, INTENT(IN) :: keep
      REAL(KIND=GP), ALLOCATABLE, TARGET :: t(:)
      INTEGER :: m
      IF (.not.ALLOCATED(a)) THEN
         CALL galloc_r1(a,n)
         RETURN
      ENDIF
      IF (SIZE(a).eq.n) RETURN
      CALL galloc_r1(t,n)
      m = MIN(SIZE(a),n)
      IF (keep.and.(m.gt.0)) THEN
         t(1:m) = a(1:m)
#if defined(GHOST_GPU)
         CALL gdev_copy(C_LOC(t),C_LOC(a),nbytes_r(m))
#endif
      ENDIF
      CALL gfree_r1(a)
      CALL MOVE_ALLOC(t,a)
      END SUBROUTINE gresize_r1

      SUBROUTINE gresize_i1(a,n,keep)
      INTEGER, ALLOCATABLE, TARGET, INTENT(INOUT) :: a(:)
      INTEGER, INTENT(IN) :: n
      LOGICAL, INTENT(IN) :: keep
      INTEGER, ALLOCATABLE, TARGET :: t(:)
      INTEGER :: m
      IF (.not.ALLOCATED(a)) THEN
         CALL galloc_i1(a,n)
         RETURN
      ENDIF
      IF (SIZE(a).eq.n) RETURN
      CALL galloc_i1(t,n)
      m = MIN(SIZE(a),n)
      IF (keep.and.(m.gt.0)) THEN
         t(1:m) = a(1:m)
#if defined(GHOST_GPU)
         CALL gdev_copy(C_LOC(t),C_LOC(a),nbytes_i(m))
#endif
      ENDIF
      CALL gfree_i1(a)
      CALL MOVE_ALLOC(t,a)
      END SUBROUTINE gresize_i1

      SUBROUTINE gresize_r2(a,n1,n2,keep)
      REAL(KIND=GP), ALLOCATABLE, TARGET, INTENT(INOUT) :: a(:,:)
      INTEGER, INTENT(IN) :: n1,n2
      LOGICAL, INTENT(IN) :: keep
      REAL(KIND=GP), ALLOCATABLE, TARGET :: t(:,:)
      INTEGER :: m
      IF (.not.ALLOCATED(a)) THEN
         CALL galloc_r2(a,n1,n2)
         RETURN
      ENDIF
      IF ((SIZE(a,1).eq.n1).and.(SIZE(a,2).eq.n2)) RETURN
      IF (SIZE(a,1).ne.n1) STOP 'gresize: the leading dimension cannot change'
      CALL galloc_r2(t,n1,n2)
      m = MIN(SIZE(a,2),n2)
      IF (keep.and.(m.gt.0)) THEN
         t(:,1:m) = a(:,1:m)
#if defined(GHOST_GPU)
         CALL gdev_copy(C_LOC(t),C_LOC(a),nbytes_r(n1*m))
#endif
      ENDIF
      CALL gfree_r2(a)
      CALL MOVE_ALLOC(t,a)
      END SUBROUTINE gresize_r2

      SUBROUTINE gresize_i2(a,n1,n2,keep)
      INTEGER, ALLOCATABLE, TARGET, INTENT(INOUT) :: a(:,:)
      INTEGER, INTENT(IN) :: n1,n2
      LOGICAL, INTENT(IN) :: keep
      INTEGER, ALLOCATABLE, TARGET :: t(:,:)
      INTEGER :: m
      IF (.not.ALLOCATED(a)) THEN
         CALL galloc_i2(a,n1,n2)
         RETURN
      ENDIF
      IF ((SIZE(a,1).eq.n1).and.(SIZE(a,2).eq.n2)) RETURN
      IF (SIZE(a,1).ne.n1) STOP 'gresize: the leading dimension cannot change'
      CALL galloc_i2(t,n1,n2)
      m = MIN(SIZE(a,2),n2)
      IF (keep.and.(m.gt.0)) THEN
         t(:,1:m) = a(:,1:m)
#if defined(GHOST_GPU)
         CALL gdev_copy(C_LOC(t),C_LOC(a),nbytes_i(n1*m))
#endif
      ENDIF
      CALL gfree_i2(a)
      CALL MOVE_ALLOC(t,a)
      END SUBROUTINE gresize_i2

      PURE FUNCTION nbytes_i(n) RESULT(nb)
      INTEGER, INTENT(IN) :: n
      INTEGER(C_SIZE_T)   :: nb
      INTEGER             :: k
      nb = INT(n,C_SIZE_T)*INT(STORAGE_SIZE(k)/8,C_SIZE_T)
      END FUNCTION nbytes_i

      PURE FUNCTION nbytes_r(n) RESULT(nb)
      INTEGER, INTENT(IN) :: n
      INTEGER(C_SIZE_T)   :: nb
      REAL(KIND=GP)       :: x
      nb = INT(n,C_SIZE_T)*INT(STORAGE_SIZE(x)/8,C_SIZE_T)
      END FUNCTION nbytes_r

      PURE FUNCTION nbytes_c(n) RESULT(nb)
      INTEGER, INTENT(IN) :: n
      INTEGER(C_SIZE_T)   :: nb
      COMPLEX(KIND=GP)    :: z
      nb = INT(n,C_SIZE_T)*INT(STORAGE_SIZE(z)/8,C_SIZE_T)
      END FUNCTION nbytes_c

      END MODULE gmem
