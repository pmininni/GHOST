!=================================================================
! FFTP3D for GPUs
! Parallel Fast Fourier Transform in 3D
!
! Performs parallel real-to-complex and complex-to-real FFTs using
! MPI and, in each task, cuFFT or hipFFT on the device resident
! arrays while the fields are being worked on the device
! (gdev_active set, see src/memmgt/gdevice_mod.f90), or the FFTW
! library on the host copies otherwise. The host path is the one of
! fftp-3: a 2D transform of the (x,y) planes, a slab exchange with
! MPI derived data types, a cache-blocked transposition and a 1D
! transform along z. The device path follows the same steps with
! the exchange done by packing the blocks on the device into
! contiguous messages that are handed to a GPU-aware MPI (or staged
! through the host copies, see gexchange in fftp_mod.F90). The
! transforms are unnormalized and the layout of the output is the
! one of fftp-3 with FFTW, so both backends are interchangeable.
!
! You should use the FFTPLANS and MPIVARS modules (see the file
! 'fftp_mod.F90') in each program that calls any of the subroutines
! in this file. Also, you must create plans for the parallel FFT
! using the 'fftp3d_create_plan' subroutine.
!
! 2007 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.uba.ar 
!
! 16 Feb 2004: Performs complex FFTs in place.
!  8 Jul 2004: itype pointers only used to store datatypes of 
!              blocks in the row and column each processor is.
!  9 Jul 2004: Transposition uses data cache blocking.
! 13 Feb 2007: Transposition uses strip mining (rreddy@psc.edu)
! 25 Aug 2009: Hybrid MPI/OpenMP support (D. Rosenberg & P. Mininni)
! 30 Aug 2009: SINGLE/DOUBLE precision (D. Rosenberg & P. Mininni)
!  3 Jan 2017: Anisotropic boxes (P. Mininni)
!  2 Sep 2026: GPU version with cuFFT/hipFFT and OpenMP offload
!
! References:
! Mininni PD, Gomez DO, Mahajan SM; Astrophys. J. 619, 1019 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Phys. Scripta T116, 123 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Adv. Sp. Res. 35, 899 (2005)
! Rosenberg DL, Mininni PD, Reddy R, Pouquet A.: Atmosph.11, 178 (2020)
!=================================================================
#include "fftw_wrappers.h"
#include "gfft_defs.h"

!=================================================================
! MODULE fftp3d_gpu
!
! Device path helpers of FFTP-GPU: execution of the device FFTs,
! transpositions, packing of the exchange blocks and the exchanges
! between tasks. The array dummies are explicit-shape with TARGET:
! inside a "target data use_device_addr" region they refer to the
! device copies of the arrays, so C_LOC yields device addresses for
! the FFT library and for MPI.
!=================================================================
      MODULE fftp3d_gpu
      IMPLICIT NONE
      CONTAINS

!*****************************************************************
      SUBROUTINE fftp3d_exec_r2c_dev(plan,in,out)
!-----------------------------------------------------------------
!
! 2D real-to-complex transform of the (x,y) planes on the device
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE fftplans
      USE gfft
      USE gdevice
      IMPLICIT NONE
      TYPE(FFTPLAN), INTENT(IN) :: plan
      REAL(KIND=GP)   , INTENT(IN) , TARGET :: in (plan%nx,plan%ny,ksta:kend)
      COMPLEX(KIND=GP), INTENT(OUT), TARGET :: out(plan%nx/2+1,plan%ny,ksta:kend)
      INTEGER(C_INT) :: ier
!$omp target data use_device_addr(in,out)
      ier = gfft_exec_r2c(plan%dplanr,C_LOC(in),C_LOC(out))
!$omp end target data
      IF (ier.ne.0) PRINT *,'fftp3d_exec_r2c_dev: error ',ier
      CALL device_sync()
      END SUBROUTINE fftp3d_exec_r2c_dev

!*****************************************************************
      SUBROUTINE fftp3d_exec_c2r_dev(plan,in,out)
!-----------------------------------------------------------------
!
! 2D complex-to-real transform of the (x,y) planes on the device
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE fftplans
      USE gfft
      USE gdevice
      IMPLICIT NONE
      TYPE(FFTPLAN), INTENT(IN) :: plan
      COMPLEX(KIND=GP), INTENT(IN) , TARGET :: in (plan%nx/2+1,plan%ny,ksta:kend)
      REAL(KIND=GP)   , INTENT(OUT), TARGET :: out(plan%nx,plan%ny,ksta:kend)
      INTEGER(C_INT) :: ier
!$omp target data use_device_addr(in,out)
      ier = gfft_exec_c2r(plan%dplanr,C_LOC(in),C_LOC(out))
!$omp end target data
      IF (ier.ne.0) PRINT *,'fftp3d_exec_c2r_dev: error ',ier
      CALL device_sync()
      END SUBROUTINE fftp3d_exec_c2r_dev

!*****************************************************************
      SUBROUTINE fftp3d_exec_c2c_dev(plan,a,dir)
!-----------------------------------------------------------------
!
! 1D complex transform along z, in place, on the device
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE fftplans
      USE gfft
      USE gdevice
      IMPLICIT NONE
      TYPE(FFTPLAN), INTENT(IN) :: plan
      COMPLEX(KIND=GP), INTENT(INOUT), TARGET :: a(plan%nz,plan%ny,ista:iend)
      INTEGER, INTENT(IN) :: dir
      INTEGER(C_INT) :: ier
!$omp target data use_device_addr(a)
      ier = gfft_exec_c2c(plan%dplanc,C_LOC(a),C_LOC(a),INT(dir,C_INT))
!$omp end target data
      IF (ier.ne.0) PRINT *,'fftp3d_exec_c2c_dev: error ',ier
      CALL device_sync()
      END SUBROUTINE fftp3d_exec_c2c_dev

!*****************************************************************
      SUBROUTINE fftp3d_tra_fwd_dev(plan,c1,out)
!-----------------------------------------------------------------
!
! Transposition c1(i,j,k) -> out(k,j,i) on the device. The innermost
! loop runs along the first index of the array that is written, so
! that the writes are coalesced (154 GB/s on a MI210 against 71 GB/s
! for the other order; OpenMP tiling through team-private arrays is
! slower with the current compilers, and the transposition cannot be
! folded into a strided batched plan of the FFT library because the
! batch index would have to run over (i,j) on input and over (j,i) on
! output).
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE fftplans
      IMPLICIT NONE
      TYPE(FFTPLAN), INTENT(IN) :: plan
      COMPLEX(KIND=GP), INTENT (IN) :: c1 (ista:iend,plan%ny,plan%nz)
      COMPLEX(KIND=GP), INTENT(OUT) :: out(plan%nz,plan%ny,ista:iend)
      INTEGER :: i,j,k
!$omp target teams distribute parallel do collapse(3)
      DO i = ista,iend
         DO j = 1,plan%ny
            DO k = 1,plan%nz
               out(k,j,i) = c1(i,j,k)
            END DO
         END DO
      END DO
      END SUBROUTINE fftp3d_tra_fwd_dev

!*****************************************************************
      SUBROUTINE fftp3d_tra_bwd_dev(plan,in,c1)
!-----------------------------------------------------------------
!
! Transposition in(k,j,i) -> c1(i,j,k) on the device (see above)
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE fftplans
      IMPLICIT NONE
      TYPE(FFTPLAN), INTENT(IN) :: plan
      COMPLEX(KIND=GP), INTENT (IN) :: in(plan%nz,plan%ny,ista:iend)
      COMPLEX(KIND=GP), INTENT(OUT) :: c1(ista:iend,plan%ny,plan%nz)
      INTEGER :: i,j,k
!$omp target teams distribute parallel do collapse(3)
      DO k = 1,plan%nz
         DO j = 1,plan%ny
            DO i = ista,iend
               c1(i,j,k) = in(k,j,i)
            END DO
         END DO
      END DO
      END SUBROUTINE fftp3d_tra_bwd_dev

!*****************************************************************
      SUBROUTINE fftp3d_pack_dev(carr,sbuf,nxh,ny,ntot,i0,lr,o)
!-----------------------------------------------------------------
!
! Packs the block of carr with kx in [i0,i0+lr-1] contiguously at
! sbuf(o+1:), on the device
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nxh,ny,ntot,i0,lr,o
      COMPLEX(KIND=GP), INTENT (IN)   :: carr(nxh,ny,ksta:kend)
      COMPLEX(KIND=GP), INTENT(INOUT) :: sbuf(ntot)
      INTEGER :: i,j,k
!$omp target teams distribute parallel do collapse(3)
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,lr
               sbuf(o+i+lr*((j-1)+ny*(k-ksta))) = carr(i0+i-1,j,k)
            END DO
         END DO
      END DO
      END SUBROUTINE fftp3d_pack_dev

!*****************************************************************
      SUBROUTINE fftp3d_unpack_dev(sbuf,carr,nxh,ny,ntot,i0,lr,o)
!-----------------------------------------------------------------
!
! Unpacks the block at sbuf(o+1:) into carr, kx in [i0,i0+lr-1]
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nxh,ny,ntot,i0,lr,o
      COMPLEX(KIND=GP), INTENT   (IN) :: sbuf(ntot)
      COMPLEX(KIND=GP), INTENT(INOUT) :: carr(nxh,ny,ksta:kend)
      INTEGER :: i,j,k
!$omp target teams distribute parallel do collapse(3)
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,lr
               carr(i0+i-1,j,k) = sbuf(o+i+lr*((j-1)+ny*(k-ksta)))
            END DO
         END DO
      END DO
      END SUBROUTINE fftp3d_unpack_dev

!*****************************************************************
      SUBROUTINE fftp3d_copy_dev(dst,src,n)
!-----------------------------------------------------------------
!
! Contiguous copy of n complex numbers on the device (the block a
! task exchanges with itself never goes through MPI)
!-----------------------------------------------------------------
      USE fprecision
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      COMPLEX(KIND=GP), INTENT(OUT) :: dst(n)
      COMPLEX(KIND=GP), INTENT(IN)  :: src(n)
      INTEGER :: i
!$omp target teams distribute parallel do
      DO i = 1,n
         dst(i) = src(i)
      END DO
      END SUBROUTINE fftp3d_copy_dev

!*****************************************************************
      SUBROUTINE fftp3d_exch_contig(sendbuf,nsend,recvbuf,nrecv,so,sc,ro,rc,comm)
!-----------------------------------------------------------------
!
! Exchange of contiguous blocks between tasks with device buffers:
! the block for task r starts at sendbuf(so(r)+1) with sc(r)
! elements, the block from task r is received at recvbuf(ro(r)+1)
! with rc(r) elements. The diagonal block is copied on the device.
! The buffers are entered in a use_device_addr region and handed to
! MPI, which must be GPU-aware.
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nsend,nrecv
      COMPLEX(KIND=GP), INTENT(IN) , TARGET :: sendbuf(nsend)
      COMPLEX(KIND=GP), INTENT(OUT), TARGET :: recvbuf(nrecv)
      INTEGER, INTENT(IN) :: so(0:nprocs-1),sc(0:nprocs-1)
      INTEGER, INTENT(IN) :: ro(0:nprocs-1),rc(0:nprocs-1)
      TYPE(MPI_Comm), INTENT(IN) :: comm

      CALL fftp3d_copy_dev(recvbuf(ro(myrank)+1),sendbuf(so(myrank)+1),sc(myrank))
!$omp target data use_device_addr(sendbuf,recvbuf)
      CALL fftp3d_exch_contig_do(sendbuf,recvbuf,so,sc,ro,rc,comm)
!$omp end target data
      END SUBROUTINE fftp3d_exch_contig

!*****************************************************************
      SUBROUTINE fftp3d_exch_contig_do(sendbuf,recvbuf,so,sc,ro,rc,comm)
!-----------------------------------------------------------------
!
! The MPI calls of fftp3d_exch_contig, receiving device addresses
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      IMPLICIT NONE
      COMPLEX(KIND=GP), INTENT(IN)  :: sendbuf(*)
      COMPLEX(KIND=GP), INTENT(OUT) :: recvbuf(*)
      INTEGER, INTENT(IN) :: so(0:nprocs-1),sc(0:nprocs-1)
      INTEGER, INTENT(IN) :: ro(0:nprocs-1),rc(0:nprocs-1)
      TYPE(MPI_Comm), INTENT(IN) :: comm
      TYPE(MPI_Request) :: ireq(2*nprocs)
      INTEGER :: irank,isendTo,igetFrom

! All receives and sends are posted before waiting, so that the
! transfers to and from the different tasks proceed at the same time
      DO irank = 1,nprocs-1
         igetFrom = myrank - irank
         IF (igetFrom.lt.0) igetFrom = igetFrom + nprocs
         CALL fftp3d_irecv(recvbuf(ro(igetFrom)+1),rc(igetFrom),igetFrom,comm,ireq(irank))
      END DO
      DO irank = 1,nprocs-1
         isendTo = myrank + irank
         IF (isendTo.ge.nprocs) isendTo = isendTo - nprocs
         CALL fftp3d_isend(sendbuf(so(isendTo)+1),sc(isendTo),isendTo,comm,ireq(nprocs-1+irank))
      END DO
      IF (nprocs.gt.1) CALL MPI_WAITALL(2*(nprocs-1),ireq(1:2*(nprocs-1)),MPI_STATUSES_IGNORE,ierr)
      END SUBROUTINE fftp3d_exch_contig_do

!*****************************************************************
      SUBROUTINE fftp3d_irecv(buf,n,src,comm,req)
!
! Receives n complex numbers starting at buf (sequence association)
      USE fprecision
      USE mpivars
      USE commtypes
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n,src
      COMPLEX(KIND=GP), INTENT(OUT) :: buf(n)
      TYPE(MPI_Comm), INTENT(IN) :: comm
      TYPE(MPI_Request), INTENT(OUT) :: req
      CALL MPI_IRECV(buf,n,GC_COMPLEX,src,1,comm,req,ierr)
      END SUBROUTINE fftp3d_irecv

!*****************************************************************
      SUBROUTINE fftp3d_isend(buf,n,dst,comm,req)
      USE fprecision
      USE mpivars
      USE commtypes
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n,dst
      COMPLEX(KIND=GP), INTENT(IN) :: buf(n)
      TYPE(MPI_Comm), INTENT(IN) :: comm
      TYPE(MPI_Request), INTENT(OUT) :: req
      CALL MPI_ISEND(buf,n,GC_COMPLEX,dst,1,comm,req,ierr)
      END SUBROUTINE fftp3d_isend

!*****************************************************************
      SUBROUTINE fftp3d_exch_types_host(sendbuf,nsend,recvbuf,nrecv,stype,rtype,comm)
!-----------------------------------------------------------------
!
! Exchange with derived data types on the host copies (staged mode)
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nsend,nrecv
      COMPLEX(KIND=GP), INTENT(IN)  :: sendbuf(nsend)
      COMPLEX(KIND=GP), INTENT(OUT) :: recvbuf(nrecv)
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
      END SUBROUTINE fftp3d_exch_types_host

!*****************************************************************
      SUBROUTINE fftp3d_exchange_fwd(plan,comm)
!-----------------------------------------------------------------
!
! Slab exchange gcarr -> gc1 of the device path
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      USE fftplans
      USE gdevice
      IMPLICIT NONE
      TYPE(FFTPLAN), INTENT(IN) :: plan
      TYPE(MPI_Comm), INTENT(IN) :: comm
      INTEGER :: r,nxh,ntot

      nxh  = plan%nx/2+1
      ntot = SIZE(gsbuf)
      IF (gexchange.eq.1) THEN
         DO r = 0,nprocs-1
            CALL fftp3d_pack_dev(gcarr,gsbuf,nxh,plan%ny,ntot,gioff(r),gilen(r),gsoff(r))
         END DO
         CALL fftp3d_exch_contig(gsbuf,ntot,gc1,SIZE(gc1),gsoff,gscnt,groff,grcnt,comm)
      ELSE
         CALL gdev_update_from(C_LOC(gcarr),INT(SIZE(gcarr),C_SIZE_T)*INT(STORAGE_SIZE(gcarr)/8,C_SIZE_T))
         CALL fftp3d_exch_types_host(gcarr,SIZE(gcarr),gc1,SIZE(gc1),plan%itype1,plan%itype2,comm)
         CALL gdev_update_to(C_LOC(gc1),INT(SIZE(gc1),C_SIZE_T)*INT(STORAGE_SIZE(gc1)/8,C_SIZE_T))
      ENDIF
      END SUBROUTINE fftp3d_exchange_fwd

!*****************************************************************
      SUBROUTINE fftp3d_exchange_bwd(plan,comm)
!-----------------------------------------------------------------
!
! Slab exchange gc1 -> gcarr of the device path
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      USE fftplans
      USE gdevice
      IMPLICIT NONE
      TYPE(FFTPLAN), INTENT(IN) :: plan
      TYPE(MPI_Comm), INTENT(IN) :: comm
      INTEGER :: r,nxh,ntot

      nxh  = plan%nx/2+1
      ntot = SIZE(gsbuf)
      IF (gexchange.eq.1) THEN
         CALL fftp3d_exch_contig(gc1,SIZE(gc1),gsbuf,ntot,groff,grcnt,gsoff,gscnt,comm)
         DO r = 0,nprocs-1
            CALL fftp3d_unpack_dev(gsbuf,gcarr,nxh,plan%ny,ntot,gioff(r),gilen(r),gsoff(r))
         END DO
      ELSE
         CALL gdev_update_from(C_LOC(gc1),INT(SIZE(gc1),C_SIZE_T)*INT(STORAGE_SIZE(gc1)/8,C_SIZE_T))
         CALL fftp3d_exch_types_host(gc1,SIZE(gc1),gcarr,SIZE(gcarr),plan%itype2,plan%itype1,comm)
         CALL gdev_update_to(C_LOC(gcarr),INT(SIZE(gcarr),C_SIZE_T)*INT(STORAGE_SIZE(gcarr)/8,C_SIZE_T))
      ENDIF
      END SUBROUTINE fftp3d_exchange_bwd

      END MODULE fftp3d_gpu

!*****************************************************************
      SUBROUTINE fftp3d_init_threads(err)
!-----------------------------------------------------------------
!
! Initializes FFTW threads (used by the host path).
!
! Parameters
!     err : if zero, the initialization failed
!-----------------------------------------------------------------

!$    USE threads
      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: err

#if defined(_OPENMP)
      CALL GFFTW_INIT_THREADS(err)
      IF (err.eq.0) PRINT *,'FFTP threads initialization failed!'
#else
      err = 1
#endif

      RETURN
      END SUBROUTINE fftp3d_init_threads

!*****************************************************************
      SUBROUTINE fftp3d_create_plan(plan,n,fftdir,flags)
!-----------------------------------------------------------------
!
! Creates the host (FFTW) and device (cuFFT/hipFFT) plans in each
! task, the shared device-resident buffers of the exchange, and the
! derived data types used by the host path.
!
! Parameters
!     plan   : contains the parallel 3D plan [OUT]
!     n      : the size of the dimensions of the input array [IN]
!     fftdir : the direction of the Fourier transform [IN]
!              FFTW_FORWARD or FFTW_REAL_TO_COMPLEX (-1)
!              FFTW_BACKWARD or FFTW_COMPLEX_TO_REAL (+1)
!     flags  : flags for the FFTW [IN]
!              FFTW_ESTIMATE (sub-optimal but faster)
!              FFTW_MEASURE (optimal but slower to create plans)
!              FFTW_PATIENT AND FFTW_EXHAUSTIVE are also available
!              for extra performance, but may take a long time to
!              create plans (specially when using OpenMP)
!-----------------------------------------------------------------

      USE mpivars
      USE fftplans
      USE gfft
      USE gdevice
!$    USE threads
      USE gtimer
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n(3)
      INTEGER, INTENT(IN) :: fftdir
      INTEGER, INTENT(IN) :: flags
      TYPE(FFTPLAN), INTENT(OUT) :: plan
      INTEGER(C_INT)   :: nn(2),inemb(2),onemb(2),n1(1),nk,nb,ier
      INTEGER          :: r,is,ie,ks,ke,nxh
      CHARACTER(len=16):: senv

! Host plans, as in fftp-3
      ALLOCATE ( plan%ccarr(n(3),n(2),ista:iend)    )
      ALLOCATE ( plan%carr(n(1)/2+1,n(2),ksta:kend) )
      ALLOCATE ( plan%rarr(n(1),n(2),ksta:kend)     )
#if defined(_OPENMP)
      CALL GFFTW_PLAN_WITH_NTHREADS(nth)
#endif
      IF (fftdir.eq.FFTW_REAL_TO_COMPLEX) THEN
      CALL GFFTW_PLAN_MANY_DFT_R2C(plan%planr,2,(/n(1),n(2)/),          &
                         kend-ksta+1,plan%rarr,                       &
                         (/n(1),n(2)*(kend-ksta+1)/),1,n(1)*n(2),     &
                         plan%carr,(/n(1)/2+1,n(2)*(kend-ksta+1)/),1, &
                         (n(1)/2+1)*n(2),flags)
      ELSE
      CALL GFFTW_PLAN_MANY_DFT_C2R(plan%planr,2,(/n(1),n(2)/),          &
                         kend-ksta+1,plan%carr,                       &
                         (/n(1)/2+1,n(2)*(kend-ksta+1)/),1,           &
                         (n(1)/2+1)*n(2),plan%rarr,                   &
                         (/n(1),n(2)*(kend-ksta+1)/),1,n(1)*n(2),flags)
      ENDIF
      CALL GFFTW_PLAN_MANY_DFT(plan%planc,1,n(3),n(2)*(iend-ista+1),    &
                         plan%ccarr,(iend-ista+1)*n(2)*n(3),1,n(3),      &
                         plan%ccarr,(iend-ista+1)*n(2)*n(3),1,n(3),      &
                         fftdir,flags)
      plan%nx = n(1)
      plan%ny = n(2)
      plan%nz = n(3)
      plan%fftdir = fftdir
      ALLOCATE( plan%itype1(0:nprocs-1) )
      ALLOCATE( plan%itype2(0:nprocs-1) )
      CALL fftp3d_create_block(n,nprocs,myrank,plan%itype1,plan%itype2)

! Device plans. cuFFT and hipFFT take the dimensions in C order, the
! last one being the fastest: the 2D plan transforms the (x,y)
! planes of rarr(nx,ny,ksta:kend) into carr(nx/2+1,ny,ksta:kend),
! and the 1D plan transforms along z, the first dimension of
! ccarr(nz,ny,ista:iend), in place.
      nxh   = n(1)/2+1
      nk    = kend-ksta+1
      nb    = n(2)*(iend-ista+1)
      nn    = [INT(n(2),C_INT),INT(n(1),C_INT)]
      inemb = [INT(n(2),C_INT),INT(n(1),C_INT)]
      onemb = [INT(n(2),C_INT),INT(nxh ,C_INT)]
      n1    = [INT(n(3),C_INT)]
      IF (fftdir.eq.FFTW_REAL_TO_COMPLEX) THEN
         ier = gfft_plan_many(plan%dplanr,2_C_INT,nn,inemb,1_C_INT,     &
               INT(n(1)*n(2),C_INT),onemb,1_C_INT,INT(nxh*n(2),C_INT),  &
               INT(GFFT_TYPE_R2C,C_INT),nk)
      ELSE
         ier = gfft_plan_many(plan%dplanr,2_C_INT,nn,onemb,1_C_INT,     &
               INT(nxh*n(2),C_INT),inemb,1_C_INT,INT(n(1)*n(2),C_INT),  &
               INT(GFFT_TYPE_C2R,C_INT),nk)
      ENDIF
      IF (ier.ne.0) THEN
         PRINT *,'fftp3d_create_plan: 2D device plan failed, error ',ier
         STOP
      ENDIF
      ier = gfft_plan_many(plan%dplanc,1_C_INT,n1,n1,1_C_INT,n1(1),     &
            n1,1_C_INT,n1(1),INT(GFFT_TYPE_C2C,C_INT),nb)
      IF (ier.ne.0) THEN
         PRINT *,'fftp3d_create_plan: 1D device plan failed, error ',ier
         STOP
      ENDIF

! Shared exchange buffers and block tables (once)
      IF (.not.gbuffers_ready) THEN
         ALLOCATE( gcarr(nxh,n(2),ksta:kend) )
         ALLOCATE( gc1(ista:iend,n(2),n(3))  )
         ALLOCATE( gsbuf(nxh*n(2)*nk)        )
         gcarr = 0.0_GP; gc1 = 0.0_GP; gsbuf = 0.0_GP
         CALL gdev_alloc(C_LOC(gcarr),INT(SIZE(gcarr),C_SIZE_T)*INT(STORAGE_SIZE(gcarr)/8,C_SIZE_T))
         CALL gdev_alloc(C_LOC(gc1)  ,INT(SIZE(gc1)  ,C_SIZE_T)*INT(STORAGE_SIZE(gc1)/8  ,C_SIZE_T))
         CALL gdev_alloc(C_LOC(gsbuf),INT(SIZE(gsbuf),C_SIZE_T)*INT(STORAGE_SIZE(gsbuf)/8,C_SIZE_T))
         ALLOCATE( gioff(0:nprocs-1), gilen(0:nprocs-1) )
         ALLOCATE( gkoff(0:nprocs-1), gklen(0:nprocs-1) )
         ALLOCATE( gsoff(0:nprocs-1), gscnt(0:nprocs-1) )
         ALLOCATE( groff(0:nprocs-1), grcnt(0:nprocs-1) )
         DO r = 0,nprocs-1
            CALL range(1,nxh ,nprocs,r,is,ie)
            CALL range(1,n(3),nprocs,r,ks,ke)
            gioff(r) = is; gilen(r) = ie-is+1
            gkoff(r) = ks; gklen(r) = ke-ks+1
            gscnt(r) = gilen(r)*n(2)*nk
            grcnt(r) = (iend-ista+1)*n(2)*gklen(r)
            groff(r) = (iend-ista+1)*n(2)*(ks-1)
            IF (r.eq.0) THEN
               gsoff(r) = 0
            ELSE
               gsoff(r) = gsoff(r-1)+gscnt(r-1)
            ENDIF
         END DO
         CALL GET_ENVIRONMENT_VARIABLE('GHOST_GPU_EXCHANGE',senv)
         IF (senv(1:6).eq.'staged') THEN
            gexchange = 2
            IF (myrank.eq.0) PRINT *,'FFTP-GPU: exchange staged through the host'
         ELSE
            gexchange = 1
         ENDIF
         gbuffers_ready = .TRUE.
      ENDIF

      CALL GTInitHandle(hcom,GT_WTIME)
      CALL GTInitHandle(hfft,GT_WTIME)
      CALL GTInitHandle(htra,GT_WTIME)
      CALL GTInitHandle(htot,GT_WTIME)

      RETURN
      END SUBROUTINE fftp3d_create_plan

!*****************************************************************
      SUBROUTINE fftp3d_destroy_plan(plan)
!-----------------------------------------------------------------
!
! Destroys the plans in each task.
!
! Parameters
!     plan : the parallel 3D plan [INOUT]
!-----------------------------------------------------------------

      USE fftplans
      USE gfft
      USE gdevice
      USE gtimer
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(INOUT) :: plan
      INTEGER(C_INT) :: ier

      CALL GFFTW_DESTROY_PLAN(plan%planr)
      CALL GFFTW_DESTROY_PLAN(plan%planc)
      ier = gfft_destroy(plan%dplanr)
      ier = gfft_destroy(plan%dplanc)
      DEALLOCATE( plan%ccarr  )
      DEALLOCATE( plan%carr   )
      DEALLOCATE( plan%rarr   )
      DEALLOCATE( plan%itype1 )
      DEALLOCATE( plan%itype2 )
      IF (gbuffers_ready) THEN
         CALL gdev_free(C_LOC(gcarr))
         CALL gdev_free(C_LOC(gc1))
         CALL gdev_free(C_LOC(gsbuf))
         DEALLOCATE( gcarr, gc1, gsbuf )
         DEALLOCATE( gioff, gilen, gkoff, gklen, gsoff, gscnt, groff, grcnt )
         gbuffers_ready = .FALSE.
      ENDIF

      CALL GTFree(hcom)
      CALL GTFree(hfft)
      CALL GTFree(htra)
      CALL GTFree(htot)

      RETURN
      END SUBROUTINE fftp3d_destroy_plan

!*****************************************************************
      SUBROUTINE fftp3d_create_block(n,nprocs,myrank,itype1,itype2)
!-----------------------------------------------------------------
!
! Defines derived data types for sending and receiving 
! blocks of the 3D matrix between processors. The data 
! types are used to transpose the matrix during the FFT
! on the host path.
!
! Parameters
!     n      : the size of the dimensions of the input array [IN]
!     nprocs : the number of processors [IN]
!     myrank : the rank of the processor [IN]
!     itype1 : contains a derived data type for sending [OUT]
!     itype2 : contains a derived data type for receiving [OUT]
!-----------------------------------------------------------------

      USE commtypes
      IMPLICIT NONE

      TYPE(MPI_Datatype), INTENT(OUT), DIMENSION(0:nprocs-1) :: itype1,itype2
      INTEGER, INTENT(IN) :: n(3),nprocs
      INTEGER, INTENT(IN) :: myrank

      INTEGER :: ista,iend
      INTEGER :: ksta,kend
      INTEGER :: irank,krank
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

      RETURN
      END SUBROUTINE fftp3d_create_block

!*****************************************************************
      SUBROUTINE fftp3d_real_to_complex(plan,in,out,comm)
!-----------------------------------------------------------------
!
! Computes the 3D real-to-complex FFT in parallel. The 
! complex output has the same structure than the output 
! of the 3D FFTW, but the output is transposed. 
!
! Parameters
!     plan : the 3D plan created with fftp3d_create_plan [IN]
!     in   : real input array [IN]
!     out  : complex output array [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE fprecision
      USE mpivars
      USE fftplans
      USE gtimer
      USE gdevice
      USE fftp3d_gpu
!$    USE threads
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(OUT), TARGET, DIMENSION(plan%nz,plan%ny,ista:iend) :: out 
      REAL(KIND=GP)   , INTENT(IN) , TARGET, DIMENSION(plan%nx,plan%ny,ksta:kend) :: in

      TYPE(MPI_Request), DIMENSION(0:nprocs-1) :: ireq1,ireq2
      TYPE(MPI_Status)                         :: istatus
      TYPE(MPI_Comm)   , INTENT(IN)            :: comm
      INTEGER :: i,j,k
      INTEGER :: ii,jj,kk
      INTEGER :: irank
      INTEGER :: isendTo,igetFrom
      INTEGER :: istrip,iproc

      CALL GTStart(htot)

      IF (gdev_active) THEN
!
! Device path: 2D FFT of the (x,y) planes on the device
!
         CALL GTStart(hfft)
         CALL fftp3d_exec_r2c_dev(plan,in,gcarr)
         CALL GTStop(hfft)
!
! Exchange between tasks, with the blocks packed on the device
!
         CALL GTStart(hcom)
         CALL fftp3d_exchange_fwd(plan,comm)
         CALL GTStop(hcom)
!
! Transposition on the device
!
         CALL GTStart(htra)
         CALL fftp3d_tra_fwd_dev(plan,gc1,out)
         CALL GTStop(htra)
!
! 1D FFT along z on the device
!
         CALL GTStart(hfft)
         CALL fftp3d_exec_c2c_dev(plan,out,GFFT_FORWARD)
         CALL GTStop(hfft)

      ELSE
!
! Host path (fftp-3): 2D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
      CALL GFFTW_EXECUTE_DFT_R2C(plan%planr,in,plan%carr)
      CALL GTStop(hfft); 
!
! Transposes the result between nodes using 
! strip mining when nstrip>1 (rreddy@psc.edu)
!
      CALL GTStart(hcom)
      do iproc = 0, nprocs-1, nstrip
         do istrip=0, nstrip-1
            irank = iproc + istrip

            isendTo = myrank + irank
            if ( isendTo .ge. nprocs ) isendTo = isendTo - nprocs

            igetFrom = myrank - irank
            if ( igetFrom .lt. 0 ) igetFrom = igetFrom + nprocs
            CALL MPI_IRECV(gc1,1,plan%itype2(igetFrom),igetFrom,     & 
                          1,comm,ireq2(irank),ierr)

            CALL MPI_ISEND(plan%carr,1,plan%itype1(isendTo),isendTo, &
                          1,comm,ireq1(irank),ierr)
         enddo

         do istrip=0, nstrip-1
            irank = iproc + istrip
            CALL MPI_WAIT(ireq1(irank),istatus,ierr)
            CALL MPI_WAIT(ireq2(irank),istatus,ierr)
         enddo
      enddo
      CALL GTStop(hcom); 
!
! Cache friendly transposition
!
      CALL GTStart(htra)
!$omp parallel do collapse(3) private (i,j,k)
      DO ii = ista,iend,csize
         DO jj = 1,plan%ny,csize
            DO kk = 1,plan%nz,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(plan%ny,jj+csize-1)
               DO k = kk,min(plan%nz,kk+csize-1)
                  out(k,j,i) = gc1(i,j,k)
               END DO
               END DO
               END DO
            END DO
         END DO
      END DO
      CALL GTStop(htra)
!
! 1D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
      CALL GFFTW_EXECUTE_DFT(plan%planc,out,out)
      CALL GTStop(hfft); 

      ENDIF

      CALL GTStop(htot); 
     
      ! Update local accumulated timers:
      ffttime = GTGetTime(hfft)
      tratime = GTGetTime(htra)
      comtime = GTGetTime(hcom)
      tottime = GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp3d_real_to_complex

!*****************************************************************
      SUBROUTINE fftp3d_complex_to_real(plan,in,out,comm)
!-----------------------------------------------------------------
!
! Computes the 3D complex-to-real FFT in parallel. The 
! complex input has the same structure than the input 
! of the 3D FFTW, but should be transposed. The real 
! output has the same order than the output of the FFTW.
! The input data is destroyed during the computation.
!
! Parameters
!     plan : the 3D plan created with fftp3d_create_plan [IN]
!     in   : complex input array [IN]
!     out  : real output array [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE fprecision
      USE mpivars
      USE commtypes
      USE fftplans
      USE gtimer
      USE gdevice
      USE fftp3d_gpu
!$    USE threads
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(INOUT), TARGET, DIMENSION(plan%nz,plan%ny,ista:iend) :: in 
      REAL(KIND=GP)   , INTENT(OUT)  , TARGET, DIMENSION(plan%nx,plan%ny,ksta:kend) :: out

      TYPE(MPI_Request), DIMENSION(0:nprocs-1) :: ireq1,ireq2
      TYPE(MPI_Status)                         :: istatus
      TYPE(MPI_Comm)   , INTENT(IN)            :: comm
      INTEGER :: i,j,k
      INTEGER :: ii,jj,kk
      INTEGER :: irank
      INTEGER :: isendTo, igetFrom
      INTEGER :: istrip,iproc

      CALL GTStart(htot)

      IF (gdev_active) THEN
!
! Device path: 1D FFT along z on the device
!
         CALL GTStart(hfft)
         CALL fftp3d_exec_c2c_dev(plan,in,GFFT_INVERSE)
         CALL GTStop(hfft)
!
! Transposition on the device
!
         CALL GTStart(htra)
         CALL fftp3d_tra_bwd_dev(plan,in,gc1)
         CALL GTStop(htra)
!
! Exchange between tasks
!
         CALL GTStart(hcom)
         CALL fftp3d_exchange_bwd(plan,comm)
         CALL GTStop(hcom)
!
! 2D FFT of the (x,y) planes on the device
!
         CALL GTStart(hfft)
         CALL fftp3d_exec_c2r_dev(plan,gcarr,out)
         CALL GTStop(hfft)

      ELSE
!
! Host path (fftp-3): 1D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
      CALL GFFTW_EXECUTE_DFT(plan%planc,in,in)
      CALL GTStop(hfft); 
!
! Cache friendly transposition
!
      CALL GTStart(htra)
!$omp parallel do collapse(3) private (i,j,k)
      DO ii = ista,iend,csize
         DO jj = 1,plan%ny,csize
            DO kk = 1,plan%nz,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(plan%ny,jj+csize-1)
               DO k = kk,min(plan%nz,kk+csize-1)
                  gc1(i,j,k) = in(k,j,i)
               END DO
               END DO
               END DO
            END DO
         END DO
      END DO
      CALL GTStop(htra)
!
! Transposes the result between nodes using 
! strip mining when nstrip>1 (rreddy@psc.edu)
!
      CALL GTStart(hcom)
      do iproc = 0, nprocs-1, nstrip
         do istrip=0, nstrip-1
            irank = iproc + istrip

            isendTo = myrank + irank
            if ( isendTo .ge. nprocs ) isendTo = isendTo - nprocs

            igetFrom = myrank - irank
            if ( igetFrom .lt. 0 ) igetFrom = igetFrom + nprocs
            CALL MPI_IRECV(plan%carr,1,plan%itype1(igetFrom),igetFrom, & 
                          1,comm,ireq2(irank),ierr)
            CALL MPI_ISEND(gc1,1,plan%itype2(isendTo),isendTo, &
                          1,comm,ireq1(irank),ierr)
         enddo

         do istrip=0, nstrip-1
            irank = iproc + istrip
            CALL MPI_WAIT(ireq1(irank),istatus,ierr)
            CALL MPI_WAIT(ireq2(irank),istatus,ierr)
         enddo
      enddo
      CALL GTStop(hcom); 
!
! 2D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
      CALL GFFTW_EXECUTE_DFT_C2R(plan%planr,plan%carr,out)
      CALL GTStop(hfft); 

      ENDIF

      CALL GTStop(htot); 

      ! Update local accumulated timers:
      ffttime = GTGetTime(hfft)
      tratime = GTGetTime(htra)
      comtime = GTGetTime(hcom)
      tottime = GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp3d_complex_to_real

!*****************************************************************
      SUBROUTINE block3d(imin,imax,jmin,jmax,kmin,ista,iend, &
                        jsta,jend,ksta,kend,ioldtype,inewtype)
!-----------------------------------------------------------------
!
! Soubroutine for defining derived data types in 3D.
!
! Parameters
!     imin : the minimum value in the first dimension [IN]
!     imax : the maximum value in the first dimension [IN]
!     jmin : the minimum value in the second dimension [IN]
!     jmax : the maximum value in the second dimension [IN]
!     kmin : the minimum value in the third dimension [IN]
!     ista : start value of the block in the first dimension [IN]
!     iend : end value of the block in the first dimension [IN]
!     jsta : start value of the block in the second dimension [IN]
!     jend : end value of the block in the second dimension [IN]
!     ksta : start value of the block in the third dimension [IN]
!     kend : end value of the block in the third dimension [IN]
!     ioldtype: data type of the elements in the block [IN]
!     inewtype: the derived data type for the block [OUT]
!-----------------------------------------------------------------

      USE commtypes
      USE fftplans
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: ista,iend
      INTEGER, INTENT(IN)  :: jsta,jend
      INTEGER, INTENT(IN)  :: ksta,kend
      INTEGER, INTENT(IN)  :: imin,imax
      INTEGER, INTENT(IN)  :: jmin,jmax,kmin
      TYPE(MPI_Datatype), INTENT(IN)  :: ioldtype
      TYPE(MPI_Datatype), INTENT(OUT) :: inewtype

      INTEGER(MPI_ADDRESS_KIND) :: ilb,isize,idist
      TYPE(MPI_Datatype)        :: itemp,itemp2
      INTEGER :: ilen,jlen,klen
      INTEGER :: ierr
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

      RETURN
      END SUBROUTINE block3d

