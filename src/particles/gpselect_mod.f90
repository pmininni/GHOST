!=================================================================
! gpselect: deterministic stream compaction for the particle
! classes, on the host and on the device.
!
! gpsel_compact(n1,n2,flag,idx,m,ioff) returns in idx(1:m) the
! indices j in [n1,n2] with flag(j) /= 0, in ascending order, with
! the offset ioff added (idx(k) = j + ioff). The particles are
! processed in blocks of GPSEL_BS: one pass counts the flagged
! entries of each block, a serial scan turns the counts into
! offsets, and a second pass writes the indices of each block at
! its offset. The result does not depend on the number of threads
! or teams, so the local order of the particles after an exchange
! is reproducible. The callers set the flags with their own
! kernels (particles leaving a slab, particles of a task in the
! global database, holes left by departed particles...).
!
! The flag and idx arrays live on the device in offload builds
! (the callers allocate them with galloc); the block counts are a
! module-level scratch array grown on demand.
!=================================================================
MODULE gpselect
      USE gmem
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE
      INTEGER, PARAMETER, PUBLIC :: GPSEL_BS = 256
      INTEGER, ALLOCATABLE, TARGET, SAVE, PRIVATE :: blkcnt_(:)
      PUBLIC :: gpsel_compact, gpsel_cleanup
      PRIVATE

      CONTAINS

      SUBROUTINE gpsel_compact(n1,n2,flag,idx,m,ioff)
      INTEGER, INTENT(IN)    :: n1,n2,ioff
      INTEGER, INTENT(IN)    :: flag(n2)
      INTEGER, INTENT(INOUT) :: idx(*)
      INTEGER, INTENT(OUT)   :: m
      INTEGER                :: nblk,b,j,jlo,jhi,c,k

      m = 0
      IF (n2.lt.n1) RETURN
      nblk = (n2-n1+GPSEL_BS)/GPSEL_BS
      IF (.not.ALLOCATED(blkcnt_)) THEN
         CALL galloc(blkcnt_,nblk+1)
      ELSE IF (SIZE(blkcnt_).lt.nblk+1) THEN
         CALL gresize(blkcnt_,nblk+1,.false.)
      ENDIF
      CALL gpsel_compact_do(n1,n2,nblk,flag,idx,ioff,blkcnt_,SIZE(blkcnt_))
      ! The counts were written on the device only if the kernels ran there
      IF ( gdev_active ) CALL gupdate_from(blkcnt_)
      m = blkcnt_(nblk+1)
      END SUBROUTINE gpsel_compact

      SUBROUTINE gpsel_compact_do(n1,n2,nblk,flag,idx,ioff,cnt,ncnt)
      INTEGER, INTENT(IN)    :: n1,n2,nblk,ioff,ncnt
      INTEGER, INTENT(IN)    :: flag(n2)
      INTEGER, INTENT(INOUT) :: idx(*)
      INTEGER, INTENT(INOUT) :: cnt(ncnt)
      INTEGER                :: b,j,jlo,jhi,c,k,s

! Pass 1: number of flagged entries in each block
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active) private(j,jlo,jhi,c)
#else
!$omp parallel do private(j,jlo,jhi,c)
#endif
      DO b = 1,nblk
         jlo = n1 + (b-1)*GPSEL_BS
         jhi = MIN(n2,jlo+GPSEL_BS-1)
         c = 0
         DO j = jlo,jhi
            IF (flag(j).ne.0) c = c+1
         END DO
         cnt(b) = c
      END DO

! Exclusive scan of the block counts; cnt(nblk+1) is the total
#if defined(GHOST_GPU)
!$omp target if(target: gdev_active)
#endif
      s = 0
      DO b = 1,nblk
         c = cnt(b)
         cnt(b) = s
         s = s+c
      END DO
      cnt(nblk+1) = s
#if defined(GHOST_GPU)
!$omp end target
#endif

! Pass 2: each block writes its indices, in order, at its offset
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active) private(j,jlo,jhi,k)
#else
!$omp parallel do private(j,jlo,jhi,k)
#endif
      DO b = 1,nblk
         jlo = n1 + (b-1)*GPSEL_BS
         jhi = MIN(n2,jlo+GPSEL_BS-1)
         k = cnt(b)
         DO j = jlo,jhi
            IF (flag(j).ne.0) THEN
               k = k+1
               idx(k) = j + ioff
            ENDIF
         END DO
      END DO
      END SUBROUTINE gpsel_compact_do

      SUBROUTINE gpsel_cleanup()
      CALL gfree(blkcnt_)
      END SUBROUTINE gpsel_cleanup

END MODULE gpselect
