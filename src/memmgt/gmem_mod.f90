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
        MODULE PROCEDURE galloc_r1, galloc_r3, galloc_c1, galloc_c3
      END INTERFACE galloc
      INTERFACE galloc_ptr
        MODULE PROCEDURE galloc_r3p
      END INTERFACE galloc_ptr
      INTERFACE gfree
        MODULE PROCEDURE gfree_r1, gfree_r3, gfree_c1, gfree_c3
      END INTERFACE gfree
      INTERFACE gfree_ptr
        MODULE PROCEDURE gfree_r3p
      END INTERFACE gfree_ptr
      INTERFACE gupdate_to
        MODULE PROCEDURE gupdate_to_r1, gupdate_to_r3, gupdate_to_c1, gupdate_to_c3
      END INTERFACE gupdate_to
      INTERFACE gupdate_from
        MODULE PROCEDURE gupdate_from_r1, gupdate_from_r3, gupdate_from_c1, gupdate_from_c3
      END INTERFACE gupdate_from
      INTERFACE gcopy
        MODULE PROCEDURE gcopy_r3, gcopy_c3
      END INTERFACE gcopy

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
