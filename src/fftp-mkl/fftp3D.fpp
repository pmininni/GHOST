!=================================================================
! FFTP3D v3
! Parallel Fast Fourier Transform in 3D
!
! Performs parallel real-to-complex and complex-to-real FFTs 
! using MPI and INTEL MKL library in each node. Can offload FFTs
! to GPUs when OpenMP >= 5.0 is available, and the code is 
! compiled with offloading options.  You should use the FFTPLANS
! and MPIVARS modules (see the file 'fftp_mod.f90') in each 
! program that calls any of the subroutines in this file. Also, 
! you must create plans for the parallel FFT using the 
! 'fftp3d_create_plan' subroutine, which creates in turn 
! derived data types for message passing using the subroutine 
! 'fftp3d_create_block'.
!
! 2022 Facundo Pugliese & Pablo D. Mininni
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!=================================================================

!*****************************************************************
      SUBROUTINE fftp3d_init_threads(err)
!-----------------------------------------------------------------
!
! Initializes FFTW threads.
!
! Parameters
!     err : if zero, the initialization failed
!-----------------------------------------------------------------

!$    USE threads
!$    USE offloading
      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: err

      err = 0
!$    IF (nth.gt.1) THEN
!$       err = 1
!$       fft_threads = .true.  
!$    ENDIF
#if defined(DO_HYBRIDoffl)
!$    err = 0
!$    fft_threads = .false.
#endif
      
      RETURN
      END SUBROUTINE fftp3d_init_threads

!*****************************************************************
      SUBROUTINE fftp3d_create_plan(plan,n,fftdir,flags)
!-----------------------------------------------------------------
!
! Creates plans for the FFTW in each node.
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
!$    USE threads
!$    USE offloading
      USE gtimer
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n(3)
      INTEGER, INTENT(IN) :: fftdir
      INTEGER, INTENT(IN) :: flags
      TYPE(FFTPLAN), INTENT(OUT) :: plan
      
      plan%planr => null()
      plan%planc => null()
      dfterr = DftiCreateDescriptor(plan%planr,GP_FFT_PREC,DFTI_REAL,   &
                       2,[n(1),n(2)])
      dfterr = DftiSetValue(plan%planr,DFTI_CONJUGATE_EVEN_STORAGE,     &
                       DFTI_COMPLEX_COMPLEX)
      dfterr = DftiSetValue(plan%planr,DFTI_NUMBER_OF_TRANSFORMS,       &
                       kend-ksta+1)
!$    IF (fft_threads) dfterr = DftiSetValue(plan%planr,DFTI_THREAD_LIMIT,nth)
      IF (fftdir.eq.FFTW_REAL_TO_COMPLEX) THEN
        dfterr = DftiSetValue(plan%planr,DFTI_INPUT_STRIDES,            &
                       [0,1,2*(n(1)/2+1)])
        dfterr = DftiSetValue(plan%planr,DFTI_INPUT_DISTANCE,           &
                       2*(n(1)/2+1)*n(2))
        dfterr = DftiSetValue(plan%planr,DFTI_OUTPUT_STRIDES,           &
                       [0,1,   n(1)/2+1 ])
        dfterr = DftiSetValue(plan%planr,DFTI_OUTPUT_DISTANCE,          &
                       (n(1)/2+1)  *n(2))
      ELSE
        dfterr = DftiSetValue(plan%planr,DFTI_INPUT_STRIDES,            &
                       [0,1,   n(1)/2+1 ])
        dfterr = DftiSetValue(plan%planr,DFTI_INPUT_DISTANCE,           &
                       (n(1)/2+1)  *n(2))
        dfterr = DftiSetValue(plan%planr,DFTI_OUTPUT_STRIDES,           &
                       [0,1,2*(n(1)/2+1)])
        dfterr = DftiSetValue(plan%planr,DFTI_OUTPUT_DISTANCE,          &
                       2*(n(1)/2+1)*n(2))
      ENDIF
#if defined(DO_HYBRIDoffl)
!$omp dispatch device(targetdev)
#endif
      dfterr = DftiCommitDescriptor(plan%planr)
      dfterr = DftiCreateDescriptor(plan%planc,GP_FFT_PREC,DFTI_COMPLEX,&
                       1,[n(3)])
      dfterr = DftiSetValue(plan%planc,DFTI_PLACEMENT,DFTI_INPLACE)
      dfterr = DftiSetValue(plan%planc,DFTI_CONJUGATE_EVEN_STORAGE,     &
                       DFTI_COMPLEX_COMPLEX)
      dfterr = DftiSetValue(plan%planc,DFTI_NUMBER_OF_TRANSFORMS,       &
                       n(2)*(iend-ista+1))
!$    IF (fft_threads) dfterr = DftiSetValue(plan%planc,DFTI_THREAD_LIMIT,nth)
      dfterr = DftiSetValue(plan%planc,DFTI_INPUT_STRIDES,              &
                       [0,1,n(3)])
      dfterr = DftiSetValue(plan%planc,DFTI_INPUT_DISTANCE,n(3))
      dfterr = DftiSetValue(plan%planc,DFTI_OUTPUT_STRIDES,             &
                       [0,1,n(3)])
      dfterr = DftiSetValue(plan%planc,DFTI_OUTPUT_DISTANCE,n(3))
#if defined(DO_HYBRIDoffl)
!$omp dispatch device(targetdev)
#endif
      dfterr = DftiCommitDescriptor(plan%planc)
     
      plan%nx = n(1)
      plan%ny = n(2)
      plan%nz = n(3)
      ALLOCATE( plan%itype1(0:nprocs-1) )
      ALLOCATE( plan%itype2(0:nprocs-1) )
      CALL fftp3d_create_block(n,nprocs,myrank,plan%itype1, &
                              plan%itype2)

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
! Destroys FFTW plans in each node.
!
! Parameters
!     plan : the parallel 3D plan [INOUT]
!-----------------------------------------------------------------

      USE fftplans
      USE gtimer
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(INOUT) :: plan

      dfterr = DftiFreeDescriptor(plan%planr)
      dfterr = DftiFreeDescriptor(plan%planc)
      DEALLOCATE( plan%itype1 )
      DEALLOCATE( plan%itype2 )

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
! types are used to transpose the matrix during the FFT.
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

      INTEGER, INTENT(OUT), DIMENSION(0:nprocs-1) :: itype1,itype2
      INTEGER, INTENT(IN) :: n(3),nprocs
      INTEGER, INTENT(IN) :: myrank

      INTEGER :: ista,iend
      INTEGER :: ksta,kend
      INTEGER :: irank,krank
      INTEGER :: itemp1,itemp2

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
! of 3D FFTW, but the output is transposed. 
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
!$    USE threads
!$    USE offloading
      USE, INTRINSIC :: iso_c_binding
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(OUT), TARGET, DIMENSION(plan%nz,plan%ny,ista:iend) :: out
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%ny,plan%nz) :: c1
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION(:,:,:)        :: carr
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:)            :: p_carr
      REAL(KIND=GP), INTENT(IN), TARGET, DIMENSION(plan%nx,plan%ny,ksta:kend)     :: in
      INTEGER, DIMENSION(0:nprocs-1)      :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)                 :: comm
      INTEGER :: i,j,k
      INTEGER :: ii,jj,kk
      INTEGER :: irank
      INTEGER :: isendTo,igetFrom
      INTEGER :: istrip,iproc

!
! 2D FFT in each node using the FFTW library
!
      CALL GTStart(htot)
      ALLOCATE( carr(plan%nx/2+1,plan%ny,ksta:kend) )
      CALL GTStart(hfft)
      CALL C_F_POINTER(C_LOC(in),p_carr,shape=[plan%nx/2,plan%ny,kend-ksta+1])
      carr(1:plan%nx/2,:,:) = p_carr
#if defined(DO_HYBRIDoffl)
!$omp target data map(to:plan%planr) map(tofrom:carr) device(targetdev)
!$omp dispatch device(targetdev)
#endif
      dfterr = DftiComputeForward(plan%planr,carr(:,1,ksta))
#if defined(DO_HYBRIDoffl)
!$omp end target data
#endif
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
            CALL MPI_IRECV(c1,1,plan%itype2(igetFrom),igetFrom,      & 
                          1,comm,ireq2(irank),ierr)

            CALL MPI_ISEND(carr,1,plan%itype1(isendTo),isendTo, &
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
#if defined(DO_HYBRIDoffl)
!$omp target data map(to:plan%planc,ista,iend,plan%ny,plan%nz) map(from:out) device(targetdev)
!$omp target teams loop collapse(3) private(i,j,k) device(targetdev)
#else
!$omp parallel loop collapse(3) private (i,j,k)
#endif
      DO ii = ista,iend,csize
         DO jj = 1,plan%ny,csize
            DO kk = 1,plan%nz,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(plan%ny,jj+csize-1)
               DO k = kk,min(plan%nz,kk+csize-1)
                  out(k,j,i) = c1(i,j,k)
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
#if defined(DO_HYBRIDoffl)
!$omp dispatch device(targetdev)
#endif
      dfterr = DftiComputeForward(plan%planc,out(:,1,ista))
      CALL GTStop(hfft);
#if defined(DO_HYBRIDoffl)
!$omp end target data
#endif
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
!$    USE threads
!$    USE offloading
      USE, INTRINSIC :: iso_c_binding
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(plan%nz,plan%ny,ista:iend) :: in
      REAL(KIND=GP), INTENT(OUT), DIMENSION(plan%nx,plan%ny,ksta:kend)      :: out
      COMPLEX(KIND=GP), ALLOCATABLE, TARGET, DIMENSION(:,:,:)               :: carr
      REAL(KIND=GP), POINTER, DIMENSION(:,:,:)                              :: rarr
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%ny,plan%nz)                :: c1

      DOUBLE PRECISION                    :: t0, t1

      INTEGER, DIMENSION(0:nprocs-1)      :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)                 :: comm
      INTEGER :: i,j,k
      INTEGER :: ii,jj,kk
      INTEGER :: irank
      INTEGER :: isendTo, igetFrom
      INTEGER :: istrip,iproc

!
! 1D FFT in each node using the FFTW library

      CALL GTStart(htot)

      ALLOCATE( carr(plan%nx/2+1,plan%ny,ksta:kend) )
      CALL GTStart(hfft)
#if defined(DO_HYBRIDoffl)
!$omp target data map(to:in,plan%planr,plan%planc,ista,iend,plan%ny,plan%nz) map(from:c1) device(targetdev)
!$omp dispatch device(targetdev)
#endif
      dfterr = DftiComputeBackward(plan%planc,in(:,1,ista))
      CALL GTStop(hfft); 
!
! Cache friendly transposition
!
      CALL GTStart(htra)

#if defined(DO_HYBRIDoffl)
!$omp target teams loop collapse(3) private(i,j,k) device(targetdev)
#else
!$omp parallel loop collapse(3) private (i,j,k)
#endif
      DO ii = ista,iend,csize
         DO jj = 1,plan%ny,csize
            DO kk = 1,plan%nz,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(plan%ny,jj+csize-1)
               DO k = kk,min(plan%nz,kk+csize-1)
                  c1(i,j,k) = in(k,j,i)
               END DO
               END DO
               END DO
            END DO
         END DO
      END DO
#if defined(DO_HYBRIDoffl)
!$omp end target data
#endif
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
            CALL MPI_IRECV(carr,1,plan%itype1(igetFrom),igetFrom, & 
                          1,comm,ireq2(irank),ierr)
            CALL MPI_ISEND(c1,1,plan%itype2(isendTo),isendTo, &
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
#if defined(DO_HYBRIDoffl)
!$omp target data map(to:plan%planr) map(tofrom:carr) device(targetdev)
!$omp dispatch device(targetdev)
#endif
      dfterr = DftiComputeBackward(plan%planr,carr(:,1,ksta))
#if defined(DO_HYBRIDoffl)
!$omp end target data
#endif
      CALL C_F_POINTER(C_LOC(carr), rarr,                     &
                       shape=[2*(plan%nx/2+1),plan%ny,kend-ksta+1])
      out = rarr(1:plan%nx,:,:)
      CALL GTStop(hfft); 
      
      DEALLOCATE( carr )

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
      INTEGER, INTENT(IN)  :: ioldtype
      INTEGER, INTENT(OUT) :: inewtype

      INTEGER(KIND=MPI_ADDRESS_KIND) :: ilb,isize,idist
      INTEGER :: ilen,jlen,klen
      INTEGER :: itemp,itemp2
      INTEGER :: ierr

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
      CALL MPI_TYPE_CREATE_STRUCT(1,1,idist,itemp2,inewtype,ierr)
      CALL MPI_TYPE_FREE(itemp2,ierr)
      CALL MPI_TYPE_COMMIT(inewtype,ierr)

      RETURN
      END SUBROUTINE block3d
