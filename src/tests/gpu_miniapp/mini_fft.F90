!=================================================================
! mini_fft.F90
! Parallel 3D FFT of the mini-app: GHOST's FFTP algorithm (local 2D
! FFT, slab exchange, local transpose, local 1D FFT) with the FFT
! library reached through ISO_C_BINDING. The host path uses FFTW
! and is the reference; the device path uses cuFFT or hipFFT on
! arrays resident on the device, handing the library the device
! addresses obtained with use_device_addr.
!
! Three implementations of the slab exchange are selectable at run
! time with exchange_mode:
!   1: MPI derived datatypes on device buffers (GPU-aware MPI)
!   2: device-side pack into contiguous buffers, contiguous MPI
!   3: host staging: update device->host, MPI on host, update back
! In a host build modes 1 and 2 run on host buffers and mode 3 is
! identical to 1.
!=================================================================
#include "gpu_defs.h"

  MODULE mini_fft
      USE ISO_C_BINDING
      USE mini_prec
      USE mini_grid
      USE mini_mpivars
      USE mini_device
      USE mini_xfer
      IMPLICIT NONE

      INTEGER, PARAMETER :: FFTW_ESTIMATE  = 64
      INTEGER, PARAMETER :: FFTW_UNALIGNED = 2
      INTEGER, PARAMETER :: FFT_FORWARD  = -1
      INTEGER, PARAMETER :: FFT_BACKWARD =  1
      INTEGER, PARAMETER :: csize = 8
      INTEGER, SAVE :: exchange_mode = 1
      DOUBLE PRECISION, SAVE :: tfft = 0.0D0, tcom = 0.0D0
      DOUBLE PRECISION, SAVE :: ttra = 0.0D0, ttot = 0.0D0
      DOUBLE PRECISION, SAVE :: tpack = 0.0D0, tself = 0.0D0, tmpi = 0.0D0

      TYPE, PUBLIC :: FFTPLAN
        TYPE(C_PTR) :: hplanr = C_NULL_PTR      ! host 2D plan (FFTW)
        TYPE(C_PTR) :: hplanc = C_NULL_PTR      ! host 1D plan (FFTW)
#if defined(GHOST_GPU)
        GFFT_HANDLE_T :: dplanr                  ! device 2D plan
        GFFT_HANDLE_T :: dplanc                  ! device 1D plan
#endif
        COMPLEX(KIND=GP), ALLOCATABLE :: carr(:,:,:)  ! (nx/2+1,ny,ksta:kend)
        COMPLEX(KIND=GP), ALLOCATABLE :: c1(:,:,:)    ! (ista:iend,ny,nz)
        COMPLEX(KIND=GP), ALLOCATABLE :: ccarr(:,:,:) ! (nz,ny,ista:iend)
        REAL(KIND=GP)   , ALLOCATABLE :: rarr(:,:,:)  ! (nx,ny,ksta:kend)
        COMPLEX(KIND=GP), ALLOCATABLE :: sbuf(:)      ! packed send/recv
        INTEGER, ALLOCATABLE :: ioff(:),ilen(:),koff(:),klen(:)
        INTEGER, ALLOCATABLE :: soff(:),scnt(:),roff(:),rcnt(:)
        TYPE(MPI_Datatype), ALLOCATABLE :: itype1(:),itype2(:)
        INTEGER :: nx,ny,nz
      END TYPE FFTPLAN

! ---------------------- FFTW (host) interfaces -------------------
      INTERFACE
        TYPE(C_PTR) FUNCTION fftw_plan_r2c(rank,n,howmany,in,inembed,   &
                    istride,idist,out,onembed,ostride,odist,flags)      &
                    BIND(C,NAME=C_FFTW_PLAN_MANY_R2C)
          IMPORT :: C_PTR,C_INT
          INTEGER(C_INT), VALUE :: rank,howmany,istride,idist
          INTEGER(C_INT), VALUE :: ostride,odist,flags
          INTEGER(C_INT), INTENT(IN) :: n(*),inembed(*),onembed(*)
          TYPE(C_PTR), VALUE :: in,out
        END FUNCTION fftw_plan_r2c
        TYPE(C_PTR) FUNCTION fftw_plan_c2r(rank,n,howmany,in,inembed,   &
                    istride,idist,out,onembed,ostride,odist,flags)      &
                    BIND(C,NAME=C_FFTW_PLAN_MANY_C2R)
          IMPORT :: C_PTR,C_INT
          INTEGER(C_INT), VALUE :: rank,howmany,istride,idist
          INTEGER(C_INT), VALUE :: ostride,odist,flags
          INTEGER(C_INT), INTENT(IN) :: n(*),inembed(*),onembed(*)
          TYPE(C_PTR), VALUE :: in,out
        END FUNCTION fftw_plan_c2r
        TYPE(C_PTR) FUNCTION fftw_plan_c2c(rank,n,howmany,in,inembed,   &
                    istride,idist,out,onembed,ostride,odist,sign,flags) &
                    BIND(C,NAME=C_FFTW_PLAN_MANY_C2C)
          IMPORT :: C_PTR,C_INT
          INTEGER(C_INT), VALUE :: rank,howmany,istride,idist
          INTEGER(C_INT), VALUE :: ostride,odist,sign,flags
          INTEGER(C_INT), INTENT(IN) :: n(*),inembed(*),onembed(*)
          TYPE(C_PTR), VALUE :: in,out
        END FUNCTION fftw_plan_c2c
        SUBROUTINE fftw_exec_r2c(p,in,out) BIND(C,NAME=C_FFTW_EXEC_R2C)
          IMPORT :: C_PTR
          TYPE(C_PTR), VALUE :: p,in,out
        END SUBROUTINE fftw_exec_r2c
        SUBROUTINE fftw_exec_c2r(p,in,out) BIND(C,NAME=C_FFTW_EXEC_C2R)
          IMPORT :: C_PTR
          TYPE(C_PTR), VALUE :: p,in,out
        END SUBROUTINE fftw_exec_c2r
        SUBROUTINE fftw_exec_c2c(p,in,out) BIND(C,NAME=C_FFTW_EXEC_C2C)
          IMPORT :: C_PTR
          TYPE(C_PTR), VALUE :: p,in,out
        END SUBROUTINE fftw_exec_c2c
        SUBROUTINE fftw_destroy(p) BIND(C,NAME=C_FFTW_DESTROY)
          IMPORT :: C_PTR
          TYPE(C_PTR), VALUE :: p
        END SUBROUTINE fftw_destroy
      END INTERFACE

! ---------------------- cuFFT / hipFFT interfaces ----------------
#if defined(GHOST_GPU)
      INTERFACE
        INTEGER(C_INT) FUNCTION gfft_plan_many(plan,rank,n,inembed,     &
                    istride,idist,onembed,ostride,odist,ftype,batch)    &
                    BIND(C,NAME=C_GFFT_PLAN_MANY)
          IMPORT :: C_PTR,C_INT
          GFFT_HANDLE_T :: plan
          INTEGER(C_INT), VALUE :: rank,istride,idist,ostride,odist
          INTEGER(C_INT), VALUE :: ftype,batch
          INTEGER(C_INT), INTENT(IN) :: n(*),inembed(*),onembed(*)
        END FUNCTION gfft_plan_many
        INTEGER(C_INT) FUNCTION gfft_exec_r2c(plan,in,out)              &
                    BIND(C,NAME=C_GFFT_EXEC_R2C)
          IMPORT :: C_PTR,C_INT
          GFFT_HANDLE_T, VALUE :: plan
          TYPE(C_PTR), VALUE :: in,out
        END FUNCTION gfft_exec_r2c
        INTEGER(C_INT) FUNCTION gfft_exec_c2r(plan,in,out)              &
                    BIND(C,NAME=C_GFFT_EXEC_C2R)
          IMPORT :: C_PTR,C_INT
          GFFT_HANDLE_T, VALUE :: plan
          TYPE(C_PTR), VALUE :: in,out
        END FUNCTION gfft_exec_c2r
        INTEGER(C_INT) FUNCTION gfft_exec_c2c(plan,in,out,dir)          &
                    BIND(C,NAME=C_GFFT_EXEC_C2C)
          IMPORT :: C_PTR,C_INT
          GFFT_HANDLE_T, VALUE :: plan
          TYPE(C_PTR), VALUE :: in,out
          INTEGER(C_INT), VALUE :: dir
        END FUNCTION gfft_exec_c2c
        INTEGER(C_INT) FUNCTION gfft_destroy(plan) BIND(C,NAME=C_GFFT_DESTROY)
          IMPORT :: C_PTR,C_INT
          GFFT_HANDLE_T, VALUE :: plan
        END FUNCTION gfft_destroy
      END INTERFACE
#endif

    CONTAINS

!*****************************************************************
      SUBROUTINE fftp3d_create_plan(plan,n,fftdir)
!
! Creates the host and device plans and the exchange buffers, and
! maps the buffers to the device. This is the third allocation
! choke point (states, workspace, FFT buffers).
      TYPE(FFTPLAN), INTENT(OUT), TARGET :: plan
      INTEGER, INTENT(IN) :: n(3),fftdir
      INTEGER(C_INT) :: nn(2),inemb(2),onemb(2),n1(1)
      INTEGER(C_INT) :: nk,nb,flags,ier
      INTEGER :: r,is,ie,ks,ke,nxh

      plan%nx = n(1); plan%ny = n(2); plan%nz = n(3)
      nxh = n(1)/2+1
      nk  = kend-ksta+1
      nb  = n(2)*(iend-ista+1)
      ALLOCATE( plan%carr (nxh,n(2),ksta:kend)  )
      ALLOCATE( plan%rarr (n(1),n(2),ksta:kend) )
      ALLOCATE( plan%c1   (ista:iend,n(2),n(3)) )
      ALLOCATE( plan%ccarr(n(3),n(2),ista:iend) )
      ALLOCATE( plan%sbuf (nxh*n(2)*nk) )
      plan%carr = 0.0_GP; plan%rarr = 0.0_GP; plan%c1 = 0.0_GP
      plan%ccarr = 0.0_GP; plan%sbuf = 0.0_GP

! Block tables: ioff/ilen are the kx ranges of every rank, koff/klen
! the z ranges. soff/scnt are the packed send blocks (to rank r),
! roff/rcnt the contiguous receive blocks in c1 (from rank r).
      ALLOCATE( plan%ioff(0:nprocs-1), plan%ilen(0:nprocs-1) )
      ALLOCATE( plan%koff(0:nprocs-1), plan%klen(0:nprocs-1) )
      ALLOCATE( plan%soff(0:nprocs-1), plan%scnt(0:nprocs-1) )
      ALLOCATE( plan%roff(0:nprocs-1), plan%rcnt(0:nprocs-1) )
      DO r = 0,nprocs-1
         CALL range(1,nxh ,nprocs,r,is,ie)
         CALL range(1,n(3),nprocs,r,ks,ke)
         plan%ioff(r) = is; plan%ilen(r) = ie-is+1
         plan%koff(r) = ks; plan%klen(r) = ke-ks+1
         plan%scnt(r) = plan%ilen(r)*n(2)*nk
         plan%rcnt(r) = (iend-ista+1)*n(2)*plan%klen(r)
         plan%roff(r) = (iend-ista+1)*n(2)*(ks-1)
         IF (r.eq.0) THEN
            plan%soff(r) = 0
         ELSE
            plan%soff(r) = plan%soff(r-1)+plan%scnt(r-1)
         ENDIF
      END DO
      CALL dev_alloc_c(plan%carr,SIZE(plan%carr))
      CALL dev_alloc_r(plan%rarr,SIZE(plan%rarr))
      CALL dev_alloc_c(plan%c1  ,SIZE(plan%c1  ))
      CALL dev_alloc_c(plan%sbuf,SIZE(plan%sbuf))

! Host plans. C order: the last dimension is the fastest one.
      flags = FFTW_ESTIMATE + FFTW_UNALIGNED
      nn    = [INT(n(2),C_INT),INT(n(1),C_INT)]
      inemb = [INT(n(2),C_INT),INT(n(1),C_INT)]
      onemb = [INT(n(2),C_INT),INT(nxh ,C_INT)]
      n1    = [INT(n(3),C_INT)]
      IF (fftdir.eq.FFT_FORWARD) THEN
         plan%hplanr = fftw_plan_r2c(2_C_INT,nn,nk,C_LOC(plan%rarr),    &
                       inemb,1_C_INT,INT(n(1)*n(2),C_INT),              &
                       C_LOC(plan%carr),onemb,1_C_INT,                  &
                       INT(nxh*n(2),C_INT),flags)
      ELSE
         plan%hplanr = fftw_plan_c2r(2_C_INT,nn,nk,C_LOC(plan%carr),    &
                       onemb,1_C_INT,INT(nxh*n(2),C_INT),               &
                       C_LOC(plan%rarr),inemb,1_C_INT,                  &
                       INT(n(1)*n(2),C_INT),flags)
      ENDIF
      plan%hplanc = fftw_plan_c2c(1_C_INT,n1,nb,C_LOC(plan%ccarr),n1,   &
                    1_C_INT,n1(1),C_LOC(plan%ccarr),n1,1_C_INT,n1(1),   &
                    INT(fftdir,C_INT),flags)
      IF (.not.C_ASSOCIATED(plan%hplanr).or..not.C_ASSOCIATED(plan%hplanc)) &
         STOP 'fftp3d_create_plan: FFTW planning failed'

! Device plans
#if defined(GHOST_GPU)
      IF (fftdir.eq.FFT_FORWARD) THEN
         ier = gfft_plan_many(plan%dplanr,2_C_INT,nn,inemb,1_C_INT,     &
               INT(n(1)*n(2),C_INT),onemb,1_C_INT,INT(nxh*n(2),C_INT),  &
               INT(GFFT_TYPE_R2C,C_INT),nk)
      ELSE
         ier = gfft_plan_many(plan%dplanr,2_C_INT,nn,onemb,1_C_INT,     &
               INT(nxh*n(2),C_INT),inemb,1_C_INT,INT(n(1)*n(2),C_INT),  &
               INT(GFFT_TYPE_C2R,C_INT),nk)
      ENDIF
      IF (ier.ne.0) PRINT *,'fftp3d_create_plan: 2D device plan failed ',ier
      ier = gfft_plan_many(plan%dplanc,1_C_INT,n1,n1,1_C_INT,n1(1),     &
            n1,1_C_INT,n1(1),INT(GFFT_TYPE_C2C,C_INT),nb)
      IF (ier.ne.0) PRINT *,'fftp3d_create_plan: 1D device plan failed ',ier
#endif

! MPI derived datatypes for the exchange
      ALLOCATE( plan%itype1(0:nprocs-1), plan%itype2(0:nprocs-1) )
      CALL fftp3d_create_block(n,nprocs,myrank,plan%itype1,plan%itype2)
      END SUBROUTINE fftp3d_create_plan

!*****************************************************************
      SUBROUTINE fftp3d_destroy_plan(plan)
      TYPE(FFTPLAN), INTENT(INOUT) :: plan
      INTEGER :: r
#if defined(GHOST_GPU)
      INTEGER(C_INT) :: ier
      ier = gfft_destroy(plan%dplanr)
      ier = gfft_destroy(plan%dplanc)
#endif
      CALL dev_free_c(plan%carr,SIZE(plan%carr))
      CALL dev_free_r(plan%rarr,SIZE(plan%rarr))
      CALL dev_free_c(plan%c1  ,SIZE(plan%c1  ))
      CALL dev_free_c(plan%sbuf,SIZE(plan%sbuf))
      CALL fftw_destroy(plan%hplanr)
      CALL fftw_destroy(plan%hplanc)
      DO r = 0,nprocs-1
         CALL MPI_TYPE_FREE(plan%itype1(r),ierr)
         CALL MPI_TYPE_FREE(plan%itype2(r),ierr)
      END DO
      DEALLOCATE( plan%carr,plan%rarr,plan%c1,plan%ccarr,plan%sbuf )
      DEALLOCATE( plan%itype1,plan%itype2 )
      DEALLOCATE( plan%ioff,plan%ilen,plan%koff,plan%klen )
      DEALLOCATE( plan%soff,plan%scnt,plan%roff,plan%rcnt )
      END SUBROUTINE fftp3d_destroy_plan

!*****************************************************************
      SUBROUTINE fftp3d_create_block(n,nprocs,myrank,itype1,itype2)
!
! Verbatim from GHOST's fftp3D.fpp
      TYPE(MPI_Datatype), INTENT(OUT), DIMENSION(0:nprocs-1) :: itype1,itype2
      INTEGER, INTENT(IN) :: n(3),nprocs,myrank
      INTEGER :: ista,iend,ksta,kend,irank,krank
      TYPE(MPI_Datatype) :: itemp1,itemp2

      CALL range(1,n(3),nprocs,myrank,ksta,kend)
      DO irank = 0,nprocs-1
         CALL range(1,n(1)/2+1,nprocs,irank,ista,iend)
         CALL block3d(1,n(1)/2+1,1,n(2),ksta,ista,iend,1,n(2), &
                     ksta,kend,GC_COMPLEX,itemp1)
         itype1(irank) = itemp1
      END DO
      CALL range(1,n(1)/2+1,nprocs,myrank,ista,iend)
      DO krank = 0,nprocs-1
         CALL range(1,n(3),nprocs,krank,ksta,kend)
         CALL block3d(ista,iend,1,n(2),1,ista,iend,1,n(2),     &
                     ksta,kend,GC_COMPLEX,itemp2)
         itype2(krank) = itemp2
      END DO
      END SUBROUTINE fftp3d_create_block

!*****************************************************************
      SUBROUTINE block3d(imin,imax,jmin,jmax,kmin,ista,iend, &
                        jsta,jend,ksta,kend,ioldtype,inewtype)
!
! Verbatim from GHOST's fftp3D.fpp
      INTEGER, INTENT(IN)  :: ista,iend,jsta,jend,ksta,kend
      INTEGER, INTENT(IN)  :: imin,imax,jmin,jmax,kmin
      TYPE(MPI_Datatype), INTENT(IN)  :: ioldtype
      TYPE(MPI_Datatype), INTENT(OUT) :: inewtype
      INTEGER(MPI_ADDRESS_KIND) :: ilb,isize,idist
      TYPE(MPI_Datatype)        :: itemp,itemp2
      INTEGER :: ilen,jlen,klen,ierr
      INTEGER :: blocklengths(1)
      INTEGER(MPI_ADDRESS_KIND) :: displs(1)
      TYPE(MPI_Datatype) :: types(1)

      CALL MPI_TYPE_GET_EXTENT(ioldtype,ilb,isize,ierr)
      ilen = iend-ista+1
      jlen = jend-jsta+1
      klen = kend-ksta+1
      CALL MPI_TYPE_VECTOR(jlen,ilen,imax-imin+1,ioldtype,itemp,ierr)
      idist = (imax-imin+1)*(jmax-jmin+1)*isize
      CALL MPI_TYPE_CREATE_HVECTOR(klen,1,idist,itemp,itemp2,ierr)
      CALL MPI_TYPE_FREE(itemp,ierr)
      idist = ((imax-imin+1)*(jmax-jmin+1)*(ksta-kmin) &
              +(imax-imin+1)*(jsta-jmin)+(ista-imin))*isize
      blocklengths(1) = 1
      displs(1) = idist
      types(1) = itemp2
      CALL MPI_TYPE_CREATE_STRUCT(1,blocklengths,displs,types,inewtype,ierr)
      CALL MPI_TYPE_FREE(itemp2,ierr)
      CALL MPI_TYPE_COMMIT(inewtype,ierr)
      END SUBROUTINE block3d

! ================================================================
! Device FFT execution helpers. The dummies are explicit-shape with
! TARGET so that C_LOC is legal; inside use_device_addr they refer
! to the device copies, so C_LOC yields device addresses.
! ================================================================
      SUBROUTINE exec_r2c_dev(plan,in,out)
      TYPE(FFTPLAN), INTENT(IN) :: plan
      REAL(KIND=GP)   , INTENT(IN) , TARGET :: in (plan%nx,plan%ny,ksta:kend)
      COMPLEX(KIND=GP), INTENT(OUT), TARGET :: out(plan%nx/2+1,plan%ny,ksta:kend)
#if defined(GHOST_GPU)
      INTEGER(C_INT) :: ier
!$omp target data use_device_addr(in,out)
      ier = gfft_exec_r2c(plan%dplanr,C_LOC(in),C_LOC(out))
!$omp end target data
      IF (ier.ne.0) PRINT *,'exec_r2c_dev: error ',ier
      CALL device_sync()
#endif
      END SUBROUTINE exec_r2c_dev

      SUBROUTINE exec_c2r_dev(plan,in,out)
      TYPE(FFTPLAN), INTENT(IN) :: plan
      COMPLEX(KIND=GP), INTENT(IN) , TARGET :: in (plan%nx/2+1,plan%ny,ksta:kend)
      REAL(KIND=GP)   , INTENT(OUT), TARGET :: out(plan%nx,plan%ny,ksta:kend)
#if defined(GHOST_GPU)
      INTEGER(C_INT) :: ier
!$omp target data use_device_addr(in,out)
      ier = gfft_exec_c2r(plan%dplanr,C_LOC(in),C_LOC(out))
!$omp end target data
      IF (ier.ne.0) PRINT *,'exec_c2r_dev: error ',ier
      CALL device_sync()
#endif
      END SUBROUTINE exec_c2r_dev

      SUBROUTINE exec_c2c_dev(plan,a,dir)
      TYPE(FFTPLAN), INTENT(IN) :: plan
      COMPLEX(KIND=GP), INTENT(INOUT), TARGET :: a(plan%nz,plan%ny,ista:iend)
      INTEGER, INTENT(IN) :: dir
#if defined(GHOST_GPU)
      INTEGER(C_INT) :: ier
!$omp target data use_device_addr(a)
      ier = gfft_exec_c2c(plan%dplanc,C_LOC(a),C_LOC(a),INT(dir,C_INT))
!$omp end target data
      IF (ier.ne.0) PRINT *,'exec_c2c_dev: error ',ier
      CALL device_sync()
#endif
      END SUBROUTINE exec_c2c_dev

! ================================================================
! Local transposes. Device: flat collapsed loop. Host: GHOST's
! cache-blocked loop.
! ================================================================
      SUBROUTINE tra_fwd(c1,out,ondev)
      COMPLEX(KIND=GP), INTENT (IN) :: c1 (ista:iend,ny,nz)
      COMPLEX(KIND=GP), INTENT(OUT) :: out(nz,ny,ista:iend)
      LOGICAL, INTENT(IN) :: ondev
      INTEGER :: i,j,k,ii,jj,kk

      IF (ondev) THEN
#if defined(GHOST_GPU)
#if GPU_LOOP_SPELLING == 2
!$omp target teams loop collapse(3)
#else
!$omp target teams distribute parallel do collapse(3)
#endif
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
                  out(k,j,i) = c1(i,j,k)
               END DO
            END DO
         END DO
         RETURN
#endif
      ENDIF
!$omp parallel do collapse(3) private (i,j,k)
      DO ii = ista,iend,csize
         DO jj = 1,ny,csize
            DO kk = 1,nz,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(ny,jj+csize-1)
               DO k = kk,min(nz,kk+csize-1)
                  out(k,j,i) = c1(i,j,k)
               END DO
               END DO
               END DO
            END DO
         END DO
      END DO
      END SUBROUTINE tra_fwd

      SUBROUTINE tra_bwd(in,c1,ondev)
      COMPLEX(KIND=GP), INTENT (IN) :: in(nz,ny,ista:iend)
      COMPLEX(KIND=GP), INTENT(OUT) :: c1(ista:iend,ny,nz)
      LOGICAL, INTENT(IN) :: ondev
      INTEGER :: i,j,k,ii,jj,kk

      IF (ondev) THEN
#if defined(GHOST_GPU)
#if GPU_LOOP_SPELLING == 2
!$omp target teams loop collapse(3)
#else
!$omp target teams distribute parallel do collapse(3)
#endif
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
                  c1(i,j,k) = in(k,j,i)
               END DO
            END DO
         END DO
         RETURN
#endif
      ENDIF
!$omp parallel do collapse(3) private (i,j,k)
      DO ii = ista,iend,csize
         DO jj = 1,ny,csize
            DO kk = 1,nz,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(ny,jj+csize-1)
               DO k = kk,min(nz,kk+csize-1)
                  c1(i,j,k) = in(k,j,i)
               END DO
               END DO
               END DO
            END DO
         END DO
      END DO
      END SUBROUTINE tra_bwd

! ================================================================
! Pack / unpack kernels for exchange mode 2. The block of carr
! with kx in [i0,i0+lr-1] is packed contiguously at sbuf(o+1:).
! ================================================================
      SUBROUTINE pack_block(carr,sbuf,nxh,ntot,i0,lr,o,ondev)
      INTEGER, INTENT(IN) :: nxh,ntot,i0,lr,o
      COMPLEX(KIND=GP), INTENT (IN) :: carr(nxh,ny,ksta:kend)
      COMPLEX(KIND=GP), INTENT(INOUT) :: sbuf(ntot)
      LOGICAL, INTENT(IN) :: ondev
      INTEGER :: i,j,k

      IF (ondev) THEN
#if defined(GHOST_GPU)
#if GPU_LOOP_SPELLING == 2
!$omp target teams loop collapse(3)
#else
!$omp target teams distribute parallel do collapse(3)
#endif
         DO k = ksta,kend
            DO j = 1,ny
               DO i = 1,lr
                  sbuf(o+i+lr*((j-1)+ny*(k-ksta))) = carr(i0+i-1,j,k)
               END DO
            END DO
         END DO
         RETURN
#endif
      ENDIF
!$omp parallel do collapse(2) private (i)
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,lr
               sbuf(o+i+lr*((j-1)+ny*(k-ksta))) = carr(i0+i-1,j,k)
            END DO
         END DO
      END DO
      END SUBROUTINE pack_block

      SUBROUTINE unpack_block(sbuf,carr,nxh,ntot,i0,lr,o,ondev)
      INTEGER, INTENT(IN) :: nxh,ntot,i0,lr,o
      COMPLEX(KIND=GP), INTENT   (IN) :: sbuf(ntot)
      COMPLEX(KIND=GP), INTENT(INOUT) :: carr(nxh,ny,ksta:kend)
      LOGICAL, INTENT(IN) :: ondev
      INTEGER :: i,j,k

      IF (ondev) THEN
#if defined(GHOST_GPU)
#if GPU_LOOP_SPELLING == 2
!$omp target teams loop collapse(3)
#else
!$omp target teams distribute parallel do collapse(3)
#endif
         DO k = ksta,kend
            DO j = 1,ny
               DO i = 1,lr
                  carr(i0+i-1,j,k) = sbuf(o+i+lr*((j-1)+ny*(k-ksta)))
               END DO
            END DO
         END DO
         RETURN
#endif
      ENDIF
!$omp parallel do collapse(2) private (i)
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,lr
               carr(i0+i-1,j,k) = sbuf(o+i+lr*((j-1)+ny*(k-ksta)))
            END DO
         END DO
      END DO
      END SUBROUTINE unpack_block

! ================================================================
! MPI exchanges. exch_types uses GHOST's derived datatypes (one
! ISEND/IRECV pair per rank, waited immediately, as with nstrip=1).
! exch_contig sends packed blocks with contiguous counts. Both
! have an on-device variant: the buffers are entered in a
! use_device_addr region and passed down, so that MPI receives
! device addresses.
! ================================================================
      SUBROUTINE exch_types(sendbuf,nsend,recvbuf,nrecv,stype,rtype,comm,ondev)
      INTEGER, INTENT(IN) :: nsend,nrecv
      COMPLEX(KIND=GP), INTENT(IN) , TARGET :: sendbuf(nsend)
      COMPLEX(KIND=GP), INTENT(OUT), TARGET :: recvbuf(nrecv)
      TYPE(MPI_Datatype), INTENT(IN) :: stype(0:nprocs-1),rtype(0:nprocs-1)
      TYPE(MPI_Comm), INTENT(IN) :: comm
      LOGICAL, INTENT(IN) :: ondev

      IF (ondev) THEN
#if defined(GHOST_GPU)
!$omp target data use_device_addr(sendbuf,recvbuf)
         CALL exch_types_do(sendbuf,recvbuf,stype,rtype,comm)
!$omp end target data
         RETURN
#endif
      ENDIF
      CALL exch_types_do(sendbuf,recvbuf,stype,rtype,comm)
      END SUBROUTINE exch_types

      SUBROUTINE exch_types_do(sendbuf,recvbuf,stype,rtype,comm)
      COMPLEX(KIND=GP), INTENT(IN)  :: sendbuf(*)
      COMPLEX(KIND=GP), INTENT(OUT) :: recvbuf(*)
      TYPE(MPI_Datatype), INTENT(IN) :: stype(0:nprocs-1),rtype(0:nprocs-1)
      TYPE(MPI_Comm), INTENT(IN) :: comm
      TYPE(MPI_Request) :: ireq1,ireq2
      TYPE(MPI_Status)  :: istatus
      INTEGER :: irank,isendTo,igetFrom

      DO irank = 0,nprocs-1
         isendTo = myrank + irank
         IF (isendTo.ge.nprocs) isendTo = isendTo - nprocs
         igetFrom = myrank - irank
         IF (igetFrom.lt.0) igetFrom = igetFrom + nprocs
         CALL MPI_IRECV(recvbuf,1,rtype(igetFrom),igetFrom,1,comm,ireq2,ierr)
         CALL MPI_ISEND(sendbuf,1,stype(isendTo) ,isendTo ,1,comm,ireq1,ierr)
         CALL MPI_WAIT(ireq1,istatus,ierr)
         CALL MPI_WAIT(ireq2,istatus,ierr)
      END DO
      END SUBROUTINE exch_types_do

      SUBROUTINE exch_contig(sendbuf,nsend,recvbuf,nrecv,soff,scnt,roff,rcnt,comm,ondev)
      INTEGER, INTENT(IN) :: nsend,nrecv
      COMPLEX(KIND=GP), INTENT(IN) , TARGET :: sendbuf(nsend)
      COMPLEX(KIND=GP), INTENT(OUT), TARGET :: recvbuf(nrecv)
      INTEGER, INTENT(IN) :: soff(0:nprocs-1),scnt(0:nprocs-1)
      INTEGER, INTENT(IN) :: roff(0:nprocs-1),rcnt(0:nprocs-1)
      TYPE(MPI_Comm), INTENT(IN) :: comm
      LOGICAL, INTENT(IN) :: ondev

      DOUBLE PRECISION :: t0

! The block a rank sends to itself never goes through MPI: it is a
! plain contiguous copy (a device kernel on the device).
      t0 = MPI_WTIME()
      CALL copy_c(recvbuf(roff(myrank)+1),sendbuf(soff(myrank)+1),scnt(myrank),ondev)
      tself = tself + (MPI_WTIME()-t0)
      t0 = MPI_WTIME()
      IF (ondev) THEN
#if defined(GHOST_GPU)
!$omp target data use_device_addr(sendbuf,recvbuf)
         CALL exch_contig_do(sendbuf,recvbuf,soff,scnt,roff,rcnt,comm)
!$omp end target data
         tmpi = tmpi + (MPI_WTIME()-t0)
         RETURN
#endif
      ENDIF
      CALL exch_contig_do(sendbuf,recvbuf,soff,scnt,roff,rcnt,comm)
      tmpi = tmpi + (MPI_WTIME()-t0)
      END SUBROUTINE exch_contig

      SUBROUTINE exch_contig_do(sendbuf,recvbuf,soff,scnt,roff,rcnt,comm)
      COMPLEX(KIND=GP), INTENT(IN)  :: sendbuf(*)
      COMPLEX(KIND=GP), INTENT(OUT) :: recvbuf(*)
      INTEGER, INTENT(IN) :: soff(0:nprocs-1),scnt(0:nprocs-1)
      INTEGER, INTENT(IN) :: roff(0:nprocs-1),rcnt(0:nprocs-1)
      TYPE(MPI_Comm), INTENT(IN) :: comm
      TYPE(MPI_Request) :: ireq1,ireq2
      TYPE(MPI_Status)  :: istatus
      INTEGER :: irank,isendTo,igetFrom

      DO irank = 1,nprocs-1
         isendTo = myrank + irank
         IF (isendTo.ge.nprocs) isendTo = isendTo - nprocs
         igetFrom = myrank - irank
         IF (igetFrom.lt.0) igetFrom = igetFrom + nprocs
         CALL irecv_c(recvbuf(roff(igetFrom)+1),rcnt(igetFrom),igetFrom,comm,ireq2)
         CALL isend_c(sendbuf(soff(isendTo)+1) ,scnt(isendTo) ,isendTo ,comm,ireq1)
         CALL MPI_WAIT(ireq1,istatus,ierr)
         CALL MPI_WAIT(ireq2,istatus,ierr)
      END DO
      END SUBROUTINE exch_contig_do

      SUBROUTINE irecv_c(buf,n,src,comm,req)
!
! Receives n complex numbers starting at buf (sequence association)
      INTEGER, INTENT(IN) :: n,src
      COMPLEX(KIND=GP), INTENT(OUT) :: buf(n)
      TYPE(MPI_Comm), INTENT(IN) :: comm
      TYPE(MPI_Request), INTENT(OUT) :: req
      CALL MPI_IRECV(buf,n,GC_COMPLEX,src,1,comm,req,ierr)
      END SUBROUTINE irecv_c

      SUBROUTINE isend_c(buf,n,dst,comm,req)
      INTEGER, INTENT(IN) :: n,dst
      COMPLEX(KIND=GP), INTENT(IN) :: buf(n)
      TYPE(MPI_Comm), INTENT(IN) :: comm
      TYPE(MPI_Request), INTENT(OUT) :: req
      CALL MPI_ISEND(buf,n,GC_COMPLEX,dst,1,comm,req,ierr)
      END SUBROUTINE isend_c

! ================================================================
! Slab exchanges: forward (carr -> c1) and backward (c1 -> carr)
! ================================================================
      SUBROUTINE exchange_fwd(plan,comm,ondev)
      TYPE(FFTPLAN), INTENT(INOUT) :: plan
      TYPE(MPI_Comm), INTENT(IN) :: comm
      LOGICAL, INTENT(IN) :: ondev
      INTEGER :: r,nxh,ntot,mode
      DOUBLE PRECISION :: t0

      nxh  = plan%nx/2+1
      ntot = SIZE(plan%sbuf)
      mode = exchange_mode
      IF (.not.ondev .and. mode.eq.3) mode = 1
      SELECT CASE (mode)
      CASE (1)
         CALL exch_types(plan%carr,SIZE(plan%carr),plan%c1,SIZE(plan%c1), &
                         plan%itype1,plan%itype2,comm,ondev)
      CASE (2)
         t0 = MPI_WTIME()
         DO r = 0,nprocs-1
            CALL pack_block(plan%carr,plan%sbuf,nxh,ntot,plan%ioff(r),   &
                            plan%ilen(r),plan%soff(r),ondev)
         END DO
         tpack = tpack + (MPI_WTIME()-t0)
         CALL exch_contig(plan%sbuf,ntot,plan%c1,SIZE(plan%c1),plan%soff,&
                          plan%scnt,plan%roff,plan%rcnt,comm,ondev)
      CASE (3)
         CALL cupd_from_n(plan%carr,SIZE(plan%carr))
         CALL exch_types(plan%carr,SIZE(plan%carr),plan%c1,SIZE(plan%c1), &
                         plan%itype1,plan%itype2,comm,.FALSE.)
         CALL cupd_to_n(plan%c1,SIZE(plan%c1))
      END SELECT
      END SUBROUTINE exchange_fwd

      SUBROUTINE exchange_bwd(plan,comm,ondev)
      TYPE(FFTPLAN), INTENT(INOUT) :: plan
      TYPE(MPI_Comm), INTENT(IN) :: comm
      LOGICAL, INTENT(IN) :: ondev
      INTEGER :: r,nxh,ntot,mode
      DOUBLE PRECISION :: t0

      nxh  = plan%nx/2+1
      ntot = SIZE(plan%sbuf)
      mode = exchange_mode
      IF (.not.ondev .and. mode.eq.3) mode = 1
      SELECT CASE (mode)
      CASE (1)
         CALL exch_types(plan%c1,SIZE(plan%c1),plan%carr,SIZE(plan%carr), &
                         plan%itype2,plan%itype1,comm,ondev)
      CASE (2)
         ! send the contiguous c1 blocks, receive packed into sbuf
         CALL exch_contig(plan%c1,SIZE(plan%c1),plan%sbuf,ntot,plan%roff, &
                          plan%rcnt,plan%soff,plan%scnt,comm,ondev)
         t0 = MPI_WTIME()
         DO r = 0,nprocs-1
            CALL unpack_block(plan%sbuf,plan%carr,nxh,ntot,plan%ioff(r), &
                              plan%ilen(r),plan%soff(r),ondev)
         END DO
         tpack = tpack + (MPI_WTIME()-t0)
      CASE (3)
         CALL cupd_from_n(plan%c1,SIZE(plan%c1))
         CALL exch_types(plan%c1,SIZE(plan%c1),plan%carr,SIZE(plan%carr), &
                         plan%itype2,plan%itype1,comm,.FALSE.)
         CALL cupd_to_n(plan%carr,SIZE(plan%carr))
      END SELECT
      END SUBROUTINE exchange_bwd

! ================================================================
! The parallel 3D transforms
! ================================================================
      SUBROUTINE fftp3d_real_to_complex(plan,in,out,comm,ondev)
      TYPE(FFTPLAN), INTENT(INOUT), TARGET :: plan
      REAL(KIND=GP)   , INTENT(IN) , TARGET :: in (plan%nx,plan%ny,ksta:kend)
      COMPLEX(KIND=GP), INTENT(OUT), TARGET :: out(plan%nz,plan%ny,ista:iend)
      TYPE(MPI_Comm), INTENT(IN) :: comm
      LOGICAL, INTENT(IN) :: ondev
      DOUBLE PRECISION :: t0,t1

      t0 = MPI_WTIME()
      IF (ondev) THEN
         CALL exec_r2c_dev(plan,in,plan%carr)
      ELSE
         CALL fftw_exec_r2c(plan%hplanr,C_LOC(in),C_LOC(plan%carr))
      ENDIF
      t1 = MPI_WTIME(); tfft = tfft + (t1-t0); t0 = t1
      CALL exchange_fwd(plan,comm,ondev)
      t1 = MPI_WTIME(); tcom = tcom + (t1-t0); t0 = t1
      CALL tra_fwd(plan%c1,out,ondev)
      t1 = MPI_WTIME(); ttra = ttra + (t1-t0); t0 = t1
      IF (ondev) THEN
         CALL exec_c2c_dev(plan,out,FFT_FORWARD)
      ELSE
         CALL fftw_exec_c2c(plan%hplanc,C_LOC(out),C_LOC(out))
      ENDIF
      t1 = MPI_WTIME(); tfft = tfft + (t1-t0)
      END SUBROUTINE fftp3d_real_to_complex

      SUBROUTINE fftp3d_complex_to_real(plan,in,out,comm,ondev)
      TYPE(FFTPLAN), INTENT(INOUT), TARGET :: plan
      COMPLEX(KIND=GP), INTENT(INOUT), TARGET :: in (plan%nz,plan%ny,ista:iend)
      REAL(KIND=GP)   , INTENT(OUT)  , TARGET :: out(plan%nx,plan%ny,ksta:kend)
      TYPE(MPI_Comm), INTENT(IN) :: comm
      LOGICAL, INTENT(IN) :: ondev
      DOUBLE PRECISION :: t0,t1

      t0 = MPI_WTIME()
      IF (ondev) THEN
         CALL exec_c2c_dev(plan,in,FFT_BACKWARD)
      ELSE
         CALL fftw_exec_c2c(plan%hplanc,C_LOC(in),C_LOC(in))
      ENDIF
      t1 = MPI_WTIME(); tfft = tfft + (t1-t0); t0 = t1
      CALL tra_bwd(in,plan%c1,ondev)
      t1 = MPI_WTIME(); ttra = ttra + (t1-t0); t0 = t1
      CALL exchange_bwd(plan,comm,ondev)
      t1 = MPI_WTIME(); tcom = tcom + (t1-t0); t0 = t1
      IF (ondev) THEN
         CALL exec_c2r_dev(plan,plan%carr,out)
      ELSE
         CALL fftw_exec_c2r(plan%hplanr,C_LOC(plan%carr),C_LOC(out))
      ENDIF
      t1 = MPI_WTIME(); tfft = tfft + (t1-t0)
      END SUBROUTINE fftp3d_complex_to_real

      SUBROUTINE cupd_from_n(a,n)
      INTEGER, INTENT(IN) :: n
      COMPLEX(KIND=GP), INTENT(INOUT) :: a(n)
#if defined(GHOST_GPU)
!$omp target update from(a)
#endif
      END SUBROUTINE cupd_from_n

      SUBROUTINE cupd_to_n(a,n)
      INTEGER, INTENT(IN) :: n
      COMPLEX(KIND=GP), INTENT(IN) :: a(n)
#if defined(GHOST_GPU)
!$omp target update to(a)
#endif
      END SUBROUTINE cupd_to_n

      SUBROUTINE fft_timers_reset()
      tfft = 0.0D0; tcom = 0.0D0; ttra = 0.0D0
      tpack = 0.0D0; tself = 0.0D0; tmpi = 0.0D0
      END SUBROUTINE fft_timers_reset

!*****************************************************************
      SUBROUTINE copy_c(dst,src,n,ondev)
!
! Contiguous copy of n complex numbers, on the device or on the host
      INTEGER, INTENT(IN) :: n
      COMPLEX(KIND=GP), INTENT(OUT) :: dst(n)
      COMPLEX(KIND=GP), INTENT(IN)  :: src(n)
      LOGICAL, INTENT(IN) :: ondev
      INTEGER :: i
      IF (ondev) THEN
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do
         DO i = 1,n
            dst(i) = src(i)
         END DO
         RETURN
#endif
      ENDIF
!$omp parallel do
      DO i = 1,n
         dst(i) = src(i)
      END DO
      END SUBROUTINE copy_c

  END MODULE mini_fft
