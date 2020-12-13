!=================================================================
! FFTP3D v3
! Parallel Fast Fourier Transform in 2D
!
! Performs parallel real-to-complex and complex-to-real FFTs 
! using MPI and the CUDA library in each node. You should use 
! the FFTPLANS and MPIVARS modules (see the file 'fftp_mod.f90') 
! in each program that CALLs any of the subroutines in this 
! file. Also, you must create plans for the parallel FFT using 
! the 'fftp3d_create_plan' subroutine, which creates in turn 
! derived data types for message passing using the subroutine 
! 'fftp3d_create_block'.
!
! 2011 Duane Rosenberg & Pablo D. Mininni
!      National Center for Atmospheric Research
!      e-mail: mininni@df.uba.ar 
!
! References:
! Mininni PD, Gomez DO, Mahajan SM; Astrophys. J. 619, 1019 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Phys. Scripta T116, 123 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Adv. Sp. Res. 35, 899 (2005)
!=================================================================
#include "fftp3D.h"
#define GGPU_TRA 

!*****************************************************************
      SUBROUTINE fftp3d_create_plan(plan,n,fftdir,flags)
!-----------------------------------------------------------------
!
! Creates plans for the CUFFT in each node.
!
! Parameters
!     plan   : contains the parallel 2D plan [OUT]
!     n      : the size of the dimensions of the input array [IN]
!     fftdir : the direction of the Fourier transform [IN]
!              FFFTCU_REAL_TO_COMPLEX (-1)
!              or FFFTCU_COMPLEX_TO_REAL (+1)
!     flags  : flags for the CUFFT [IN]
!              Curently unused.
!-----------------------------------------------------------------

      USE mpivars
      USE fftplans
      USE fprecision
      USE, INTRINSIC :: iso_c_binding
      USE cuda_bindings
      USE cutypes
      USE gtimer
 
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n(3)
      INTEGER, INTENT(IN) :: fftdir
      INTEGER, INTENT(IN) :: flags
      TYPE(FFTPLAN), INTENT(OUT) :: plan

      INTEGER             :: istr, idist, ostr, odist, nrank
      INTEGER             :: i, first, last
      INTEGER             :: iret
      INTEGER             :: na(2)
      INTEGER             :: pinembed(2), ponembed(2)
      INTEGER             :: iplan

! Allocate pinned memory in the host
      plan%szccd_= max(2* n(3)     *n(2)*(iend-ista+1)*GFLOATBYTESZ,GFLOATBYTESZ)
      plan%szcd_ = max(2*(n(1)/2+1)*n(2)*(kend-ksta+1)*GFLOATBYTESZ,GFLOATBYTESZ)
      plan%szrd_ = max(   n(1)     *n(2)*(kend-ksta+1)*GFLOATBYTESZ,GFLOATBYTESZ)

      iret = cudaHostAlloc ( plan%pccarr_, plan%szccd_, cudaHostAllocPortable) 
      cudaErrChk()
      iret = cudaHostAlloc ( plan%pcarr_ , plan%szcd_ , cudaHostAllocPortable) 
      cudaErrChk()
      iret = cudaHostAlloc ( plan%prarr_ , plan%szrd_ , cudaHostAllocPortable)
      cudaErrChk()
      iret = cudaHostAlloc ( plan%pccarrt_ , plan%szccd_ , cudaHostAllocPortable)
      cudaErrChk()

      CALL c_f_pointer(plan%pccarr_ ,plan%ccarr ,(/n(3)       ,n(2),iend-ista+1/))
      CALL c_f_pointer(plan%pccarrt_,plan%ccarrt,(/iend-ista+1,n(2),n(3)       /))
      CALL c_f_pointer(plan%pcarr_  ,plan%carr  ,(/n(1)/2+1   ,n(2),kend-ksta+1/))
      CALL c_f_pointer(plan%prarr_  ,plan%rarr  ,(/n(1)       ,n(2),kend-ksta+1/))

! Allocate memory in the device
      iret = cudaMalloc(plan%cu_ccd_ , plan%szccd_)
      cudaErrChk()
      iret = cudaMalloc(plan%cu_ccd1_, plan%szccd_)
      iret = cudaMalloc(plan%cu_cd_  , plan%szcd_ )

      cudaErrChk()
      iret = cudaMalloc(plan%cu_rd_  , plan%szrd_ )
      cudaErrChk()

!
! Create streams in the GPU
      IF (streams_created.ne.1) THEN
         streams_created = 1
         DO i = 1,nstreams
            iret = cudaStreamCreate(pstream_(i))
            cudaErrChk()
         END DO
      END IF

! Create plans and auxiliary arrays for each stream
      DO i = 1,nstreams      

         CALL range(ista,iend,nstreams,i-1,first,last)
             issta (i) = first
             issnd (i) = last
         CALL range(ksta,kend,nstreams,i-1,first,last)
             kssta (i) = first
             kssnd (i) = last
         plan%str_szccd_(i) = max(2* n(3)     *n(2)*(issnd(i)-issta(i)+1) &
                                              *GFLOATBYTESZ,GFLOATBYTESZ)
         plan%str_szcd_ (i) = max(2*(n(1)/2+1)*n(2)*(kssnd(i)-kssta(i)+1) &
                                              *GFLOATBYTESZ,GFLOATBYTESZ)
         plan%str_szrd_ (i) = max(   n(1)     *n(2)*(kssnd(i)-kssta(i)+1) &
                                              *GFLOATBYTESZ,GFLOATBYTESZ)

         IF (fftdir.eq.FFTCU_REAL_TO_COMPLEX) THEN
            ! NOTE: in the following, the first dimension is the outermost
            !       and the 2nd is the innermost (contiguous) dimension:
            nrank= 2
!!          na      (1) = n(1)         ; na      (2) = n(2)              ;
!!          pinembed(1) = n(1)         ; pinembed(2) = n(2)*(kssnd(i)-kssta(i)+1);
!!          ponembed(1) = n(1)/2+1     ; ponembed(2) = n(2)*(kssnd(i)-kssta(i)+1);

            na      (2) = n(1)         ; na      (1) = n(2);
            pinembed(2) = n(1)         ; pinembed(1) = n(2)*(kssnd(i)-kssta(i)+1);
            ponembed(2) = n(1)/2+1     ; ponembed(1) = n(2)*(kssnd(i)-kssta(i)+1);

            istr        = 1            ; idist       = n(1)*n(2)         ;
            ostr        = 1            ; odist       = n(2)*(n(1)/2+1)   ;
            iret = cufftPlanMany(plan%icuplanr_(i),nrank,na,pinembed,istr,idist, &
                            ponembed,ostr,odist,GCUFFTDEFR2C,kssnd(i)-kssta(i)+1);
            cudaErrChk()
         ELSE
            nrank= 2;
!!          na      (1) = n(1)         ; na      (2) = n(2)              ;
!!          pinembed(1) = n(1)/2+1     ; pinembed(2) = n(2)*(kssnd(i)-kssta(i)+1);
!!          ponembed(1) = n(1)         ; ponembed(2) = n(2)*(kssnd(i)-kssta(i)+1);

            na      (2) = n(1)         ; na      (1) = n(2);
            pinembed(2) = n(1)/2+1     ; pinembed(1) = n(2)*(kssnd(i)-kssta(i)+1);
            ponembed(2) = n(1)         ; ponembed(1) = n(2)*(kssnd(i)-kssta(i)+1);

            istr        = 1            ; idist       = n(2)*(n(1)/2+1)   ;
            ostr        = 1            ; odist       = n(1)*n(2)         ; 
            iret = cufftPlanMany(plan%icuplanr_(i),nrank,na,pinembed,istr,idist, &
                            ponembed,ostr,odist,GCUFFTDEFC2R,kssnd(i)-kssta(i)+1);
            cudaErrChk()
         ENDIF
         nrank       = 1
         na      (1) = n(3)                         ;
         pinembed(1) = max(n(3)*n(2)*(issnd(i)-issta(i)+1),1);
         ponembed(1) = max(n(3)*n(2)*(issnd(i)-issta(i)+1),1);

         istr        = 1            ; idist      = n(3);
         ostr        = 1            ; odist      = n(3);

         iret = cufftPlanMany(plan%icuplanc_(i),nrank,na,pinembed,istr,idist,    &
              ponembed,ostr,odist,GCUFFTDEFC2C,max(n(2)*(issnd(i)-issta(i)+1),1));
         cudaErrChk()
#if 0
         IF ( iret.ne.CUFFT_SUCCESS) THEN
            WRITE(*,*)': fftp3d_create_plan: cufftPlanMany::icuplanc::c2c failed: iret=',&
                 iret,' myrank=',myrank,                                         &
                 ' na=',na(1),' pinembed=',pinembed(1),' ponembed=',ponembed(1), &
                 ' istr=',istr,' idist=',idist,' ostr=',ostr,' odist=',odist   , &
                 ' ista=',ista,' iend=',iend
           STOP 
         ENDIF
#endif

      END DO
      
      plan%nx = n(1)
      plan%ny = n(2)
      plan%nz = n(3)
      ALLOCATE( plan%itype1(0:nprocs-1) )
      ALLOCATE( plan%itype2(0:nprocs-1) )
      CALL fftp3d_create_block(n,nprocs,myrank,plan%itype1, &
                              plan%itype2)


      ! Initialize timer handles:
      CALL GTInitHandle(hcom,GT_WTIME)
      CALL GTInitHandle(hfft,GT_WTIME)
      CALL GTInitHandle(hmem,GT_WTIME)
      CALL GTInitHandle(htra,GT_WTIME)
      CALL GTInitHandle(hass,GT_WTIME)
      CALL GTInitHandle(htot,GT_WTIME)



      RETURN
      END SUBROUTINE fftp3d_create_plan

!*****************************************************************
      SUBROUTINE fftp3d_destroy_plan(plan)
!-----------------------------------------------------------------
!
! Destroys CUFFT plans in each node.
!
! Parameters
!     plan : the parallel 2D plan [INOUT]
!-----------------------------------------------------------------

      USE fftplans
      USE cuda_bindings
      USE cutypes
      USE threads
      USE gtimer

      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(INOUT) :: plan
      INTEGER                      :: iret
      INTEGER                      :: i

      IF (streams_created.eq.1) THEN
         streams_created = 0
         DO i = 1,nstreams
            iret = cudaStreamDestroy(pstream_(i))
         END DO
      ENDIF
      DO i = 1,nstreams
         iret = cufftDestroy(plan%icuplanr_(i))
         iret = cufftDestroy(plan%icuplanc_(i))
      END DO
      iret = cudaFreeHost (plan%pccarr_)
      iret = cudaFreeHost (plan%pcarr_)
      iret = cudaFreeHost (plan%prarr_)
      iret = cudaFreeHost (plan%pccarrt_)
      iret = cudaFree(plan%cu_ccd_)
      iret = cudaFree(plan%cu_ccd1_)
      iret = cudaFree(plan%cu_cd_)
      iret = cudaFree(plan%cu_rd_)
      DEALLOCATE( plan%itype1 )
      DEALLOCATE( plan%itype2 )
!
! Free time counters
      CALL GTFree(hcom); 
      CALL GTFree(hfft); 
      CALL GTFree(hmem); 
      CALL GTFree(htra); 
      CALL GTFree(hass); 
      CALL GTFree(htot);

      RETURN
      END SUBROUTINE fftp3d_destroy_plan

!*****************************************************************
      SUBROUTINE fftp3d_create_block(n,nprocs,myrank,itype1,itype2)
!-----------------------------------------------------------------
!
! Defines derived data types for sending and receiving 
! blocks of the 2D matrix between processors. The data 
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
      USE cuda_bindings
      USE cutypes
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
! of the 3D CUFFT, but the output is transposed.
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
      USE, INTRINSIC :: iso_c_binding
      USE cuda_bindings
      USE cutypes
      USE threads
      USE gtimer
      IMPLICIT NONE

      TYPE(FFTPLAN)   , INTENT (IN)                                     :: plan
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(plan%nz,plan%ny,ista:iend) :: out 
      REAL(KIND=GP), INTENT(IN), DIMENSION(plan%nx,plan%ny,ksta:kend)   :: in

      INTEGER  (C_SIZE_T)                 :: byteoffset1,  byteoffset2
      INTEGER, DIMENSION(0:nprocs-1)      :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)                 :: comm
      INTEGER :: i,iret,j,k
      INTEGER :: ii,jj,kk
      INTEGER :: irank
      INTEGER :: isendTo,igetFrom
      INTEGER :: istrip,iproc
      REAL    :: etime

!
!     NOTE: timers, hcom,hfft,hmem,htra,htot initialized in fftp3d_create_plan
      CALL GTStart(htot)
!
! 2D real-to-complex FFT in each device using the CUFFT library
      DO i = 1,nstreams ! Set streams for each FFT plan
         iret = cufftSetStream(plan%icuplanr_(i),pstream_(i));
      END DO

      CALL GTStart(hass)
!$omp parallel do if (kend-ksta.ge.nth) private (i,j)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,plan%ny
            DO i = 1,plan%nx
              plan%rarr(i,j,k-ksta+1) = in(i,j,k)
            END DO
         END DO
      END DO
      CALL GTStop(hass)

!
! Data sent to cuFFT must reside on device:
      CALL GTStart(hmem)
      DO i = 1,nstreams
         byteoffset1 = plan%nx*plan%ny*(kssta(i)-ksta)*GFLOATBYTESZ
         byteoffset2 = plan%nx*plan%ny*(kssta(i)-ksta)*GFLOATBYTESZ
         iret = cudaMemcpyAsyncOffHost2Dev(  plan%cu_rd_, & ! Dev
                                             byteoffset1, & ! OFFSET Dev
                                             plan%prarr_, & ! Host
                                             byteoffset2, & ! OFFSET Host
                           plan%str_szrd_(i), pstream_(i) )
         cudaErrChk()
      END DO
      CALL GTStop(hmem); 
 
      CALL GTStart(hfft)
      DO i = 1,nstreams
         byteoffset1 = plan%nx*plan%ny*(kssta(i)-ksta)        *GFLOATBYTESZ
         byteoffset2 = 2*(plan%nx/2+1)*plan%ny*(kssta(i)-ksta)*GFLOATBYTESZ
         iret = GCUFFTEXECOFFR2C(plan%icuplanr_(i), plan%cu_rd_, & ! Dev
                                                    byteoffset1, & ! OFFSET
                                                    plan%cu_cd_, & ! Dev
                                                    byteoffset2)   ! OFFSET
         cudaErrChk()
      END DO
      CALL GTStop(hfft); 

      CALL GTStart(hmem)
      DO i = 1,nstreams
         byteoffset1 = 2*(plan%nx/2+1)*plan%ny*(kssta(i)-ksta)*GFLOATBYTESZ
         byteoffset2 = 2*(plan%nx/2+1)*plan%ny*(kssta(i)-ksta)*GFLOATBYTESZ
         iret = cudaMemCpyAsyncOffDev2Host(  plan%pcarr_, & ! Host
                                             byteoffset1, & ! OFFSET Host
                                             plan%cu_cd_, & ! Dev
                                             byteoffset2, & ! OFFSET Dev
                         plan%str_szcd_(i), pstream_(i) )
         cudaErrChk()
      END DO
      DO i = 1,nstreams
         iret = cudaStreamSynchronize(pstream_(i))
      END DO
      CALL GTStop(hmem); 

! NOTE: If nrocs = 1, then we can carry out the transpose directly
!       on the CUDA device.
!
! Transposes the result between nodes using 
! strip mining when nstrip>1 (rreddy@psc.edu)
!
      CALL GTStart(hcom)
      DO iproc = 0, nprocs-1, nstrip
         DO istrip=0, nstrip-1
            irank = iproc + istrip

            isendTo = myrank + irank
            IF ( isendTo .GE. nprocs ) isendTo = isendTo - nprocs

            igetFrom = myrank - irank
            IF ( igetFrom .LT. 0 ) igetFrom = igetFrom + nprocs

            CALL MPI_IRECV(plan%ccarrt,1,plan%itype2(igetFrom),igetFrom,      & 
                          1,comm,ireq2(irank),ierr)
            CALL MPI_ISEND(plan%carr,1,plan%itype1(isendTo),isendTo, &
                          1,comm,ireq1(irank),ierr)
         ENDDO

         DO istrip=0, nstrip-1
            irank = iproc + istrip
            CALL MPI_WAIT(ireq1(irank),istatus,ierr)
            CALL MPI_WAIT(ireq2(irank),istatus,ierr)
         ENDDO
      ENDDO
      CALL GTStop(hcom); 

#if defined(GGPU_TRA)
      CALL GTStart(hmem)
      iret = cudaMemCpyHost2Dev(plan%cu_ccd1_, plan%pccarrt_, plan%szccd_ )
      cudaErrChk()
      CALL GTStop(hmem); 
      CALL GTStart(htra)
      CALL cuTranspose3C(plan%cu_ccd_,plan%cu_ccd1_, (iend-ista+1), &
                         plan%ny, plan%nz)

      DO i = 1,nstreams ! Set streams for each FFT plan
         iret = cufftSetStream(plan%icuplanc_(i),pstream_(i));
      END DO
      CALL GTStop(htra); 

#else

      CALL GTStart(htra)
!$omp parallel do if ((iend-ista)/csize.ge.nth) private (jj,kk,i,j,k)
       DO ii = ista,iend,csize
!$omp parallel do if ((iend-ista)/csize.lt.nth) private (kk,i,j,k)
          DO jj = 1,plan%ny,csize
             DO kk = 1,plan%nz,csize
                DO i = ii,min(iend   ,ii+csize-1)
                DO j = jj,min(plan%ny,jj+csize-1)
                DO k = kk,min(plan%nz,kk+csize-1)
                   ! Recall that ccarr is dimensioned (nz,ny,iend-ista+1),
                   ! starting at (1,1,1):
                   plan%ccarr(k,j,i-ista+1) = plan%ccarrt(i-ista+1,j,k)
                END DO
                END DO
                END DO
             END DO
          END DO
       END DO
      CALL GTStop(htra); 
!
! 1D FFT in each node using the CUFFT library
!
      CALL GTStart(hmem)
      DO i = 1,nstreams ! Set streams for each FFT plan
         iret = cufftSetStream(plan%icuplanc_(i),pstream_(i));
      END DO

      DO i = 1,nstreams
         byteoffset1 = 2*plan%nz*plan%ny*(issta(i)-ista)*GFLOATBYTESZ
         byteoffset2 = 2*plan%nz*plan%ny*(issta(i)-ista)*GFLOATBYTESZ
         iret = cudaMemcpyAsyncOffHost2Dev(  plan%cu_ccd_, & ! Dev
                                             byteoffset1 , & ! OFFSET Dev
                                             plan%pccarr_, & ! Host
                                             byteoffset2 , & ! OFFSET Host
                         plan%str_szccd_(i), pstream_(i) )
         cudaErrChk()
      END DO
      CALL GTStop(hmem); 
#endif

      CALL GTStart(hfft)
      DO i = 1,nstreams
         byteoffset1 = 2*plan%nz*plan%ny*(issta(i)-ista)*GFLOATBYTESZ
         byteoffset2 = 2*plan%nz*plan%ny*(issta(i)-ista)*GFLOATBYTESZ
         iret = GCUFFTEXECOFFC2C(plan%icuplanc_(i), plan%cu_ccd_, & ! Dev
                                                    byteoffset1 , & ! OFFSET
                                                    plan%cu_ccd_, & ! Dev
                                                    byteoffset2 , & ! OFFSET
                                           FFTCU_REAL_TO_COMPLEX)
         cudaErrChk()
      END DO
      CALL GTStop(hfft); 

      CALL GTStart(hmem)
      DO i = 1,nstreams
         byteoffset1 = 2*plan%nz*plan%ny*(issta(i)-ista)*GFLOATBYTESZ
         byteoffset2 = 2*plan%nz*plan%ny*(issta(i)-ista)*GFLOATBYTESZ
         iret = cudaMemCpyAsyncOffDev2Host(  plan%pccarr_, & ! Host
                                             byteoffset1, & ! OFFSET Host
                                             plan%cu_ccd_, & ! Dev
                                             byteoffset2, & ! OFFSET Dev
                         plan%str_szccd_(i), pstream_(i) )
         cudaErrChk()
      END DO
      DO i = 1,nstreams
         iret = cudaStreamSynchronize(pstream_(i))
      END DO
      CALL GTStop(hmem); 


      CALL GTStart(hass)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,plan%ny
            DO k = 1,plan%nz
               out(k,j,i) = plan%ccarr(k,j,i-ista+1)
            END DO
         END DO
      END DO
      CALL GTStop(hass)

      CALL GTStop(htot); 

      ! Update local accumulated timers:
      comtime = GTGetTime(hcom)
      ffttime = GTGetTime(hfft)
      memtime = GTGetTime(hmem)
      tratime = GTGetTime(htra)
      asstime = GTGetTime(hass)
      tottime = GTGetTime(htot)


      RETURN
      END SUBROUTINE fftp3d_real_to_complex


!*****************************************************************
      SUBROUTINE fftp3d_complex_to_real(plan,in,out,comm)
!-----------------------------------------------------------------
!
! Computes the 2D complex-to-real FFT in parallel. The 
! complex input has the same structure than the input 
! of the 2D CUFFT, but should be transposed. The real 
! output has the same order than the output of the CUFFT.
! The input data is destroyed during the computation.
!
! Parameters
!     plan : the 2D plan created with fftp3d_create_plan [IN]
!     in   : complex input array [IN]
!     out  : real output array [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE fprecision
      USE mpivars
      USE commtypes
      USE fftplans
      USE, INTRINSIC :: iso_c_binding
      USE cuda_bindings
      USE cutypes
      USE threads
      USE gtimer
      IMPLICIT NONE

      TYPE(FFTPLAN)   , INTENT(IN)                                     :: plan
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(plan%nz,plan%ny,ista:iend) :: in 
      REAL(KIND=GP), INTENT(OUT), DIMENSION(plan%nx,plan%ny,ksta:kend) :: out

      INTEGER  (C_SIZE_T)                 :: byteoffset1,  byteoffset2
      INTEGER, DIMENSION(0:nprocs-1)      :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)                 :: comm
      INTEGER :: i,iret,j,k
      INTEGER :: ii,jj,kk
      INTEGER :: irank
      INTEGER :: isendTo, igetFrom
      INTEGER :: istrip,iproc

!
!     NOTE: timers, hcom,hfft,hmem,htra,htot initialized in fftp3d_create_plan
      CALL GTStart(htot)
!
      DO i = 1,nstreams ! Set streams for each FFT plan
         iret = cufftSetStream(plan%icuplanc_(i), pstream_(i));
      END DO

      CALL GTStart(hass)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,plan%ny
            DO k = 1,plan%nz
               plan%ccarr(k,j,i-ista+1) = in(k,j,i)
            END DO
         END DO
      END DO
      CALL GTStop(hass)
!
! Data sent to cuFFT must reside on device:
      CALL GTStart(hmem);
      DO i = 1,nstreams
         byteoffset1 = 2*plan%nz*plan%ny*(issta(i)-ista)*GFLOATBYTESZ
         byteoffset2 = 2*plan%nz*plan%ny*(issta(i)-ista)*GFLOATBYTESZ
         iret = cudaMemcpyAsyncOffHost2Dev(  plan%cu_ccd_, & ! Dev
                                             byteoffset1 , & ! OFFSET Dev
                                             plan%pccarr_, & ! Host
                                             byteoffset2 , & ! OFFSET Host
                         plan%str_szccd_(i), pstream_(i) )
         cudaErrChk()
      END DO
      CALL GTStop(hmem); 

      CALL GTSTart(hfft)
      DO i = 1,nstreams
         byteoffset1 = 2*plan%nz*plan%ny*(issta(i)-ista)*GFLOATBYTESZ
         byteoffset2 = 2*plan%nz*plan%ny*(issta(i)-ista)*GFLOATBYTESZ
         iret = GCUFFTEXECOFFC2C(plan%icuplanc_(i), plan%cu_ccd_, & ! Dev
                                                    byteoffset1 , & ! OFFSET
                                                    plan%cu_ccd_, & ! Dev
                                                    byteoffset2 , & ! OFFSET
                                           FFTCU_COMPLEX_TO_REAL)
         cudaErrChk()
      END DO
      CALL GTStop(hfft); 


#if defined(GGPU_TRA)
      CALL GTStart(hfft); 
      DO i = 1,nstreams
         iret = cudaStreamSynchronize(pstream_(i))
      END DO
      CALL GTStop(hfft); 

      CALL GTStart(htra)
      CALL cuTranspose3C(plan%cu_ccd1_,plan%cu_ccd_,plan%nz, &
                         plan%ny,iend-ista+1)
      CALL GTStop(htra); 

      CALL GTStart(hmem);
      iret = cudaMemCpyDev2Host(plan%pccarrt_, plan%cu_ccd1_, &
                                plan%szccd_ )
      cudaErrChk()
      CALL GTStop(hmem); 
#else

      CALL GTStart(hmem);
      DO i = 1,nstreams
         byteoffset1 = 2*plan%nz*plan%ny*(issta(i)-ista)*GFLOATBYTESZ
         byteoffset2 = 2*plan%nz*plan%ny*(issta(i)-ista)*GFLOATBYTESZ
         iret = cudaMemCpyAsyncOffDev2Host(  plan%pccarr_, & ! Host
                                             byteoffset1 , & ! OFFSET Host
                                             plan%cu_ccd_, & ! Dev
                                             byteoffset2 , & ! OFFSET Dev
                           plan%str_szccd_(i), pstream_(i) )
         cudaErrChk()
      END DO
      DO i = 1,nstreams
         iret = cudaStreamSynchronize(pstream_(i))
      END DO
      CALL GTStop(hmem); 
!
! Cache friendly transposition
!
      CALL GTStart(htra);
!$omp parallel do if ((iend-ista)/csize.ge.nth) private (jj,kk,i,j,k)
      DO ii = ista,iend,csize
!$omp parallel do if ((iend-ista)/csize.lt.nth) private (kk,i,j,k)
         DO jj = 1,plan%ny,csize
            DO kk = 1,plan%nz,csize
               DO i = ii,min(iend   ,ii+csize-1)
               DO j = jj,min(plan%ny,jj+csize-1)
               DO k = kk,min(plan%nz,kk+csize-1)
                  ! Recall that ccarr is dimensioned (nz,ny,ista:iend),
                  ! starting at (1,1,1):
                  plan%ccarrt(i-ista+1,j,k) = plan%ccarr(k,j,i-ista+1)
               END DO
               END DO
               END DO
            END DO
         END DO
      END DO
      CALL GTStop(htra); 
#endif

!
! Transposes the result between nodes using 
! strip mining when nstrip>1 (rreddy@psc.edu)
!
      CALL GTStart(hcom);
      do iproc = 0, nprocs-1, nstrip
         do istrip=0, nstrip-1
            irank = iproc + istrip

            isendTo = myrank + irank
            if ( isendTo .ge. nprocs ) isendTo = isendTo - nprocs

            igetFrom = myrank - irank
            if ( igetFrom .lt. 0 ) igetFrom = igetFrom + nprocs

            CALL MPI_IRECV(plan%carr,1,plan%itype1(igetFrom),igetFrom, & 
                          1,comm,ireq2(irank),ierr)
            CALL MPI_ISEND(plan%ccarrt,1,plan%itype2(isendTo),isendTo, &
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
! 2D FFT in each node using the CUFFT library
!
      CALL GTStart(hmem);
      DO i = 1,nstreams ! Set streams for each FFT plan
         iret = cufftSetStream(plan%icuplanr_(i), pstream_(i));
      END DO
      DO i = 1,nstreams
         byteoffset1 = 2*(plan%nx/2+1)*plan%ny*(kssta(i)-ksta)*GFLOATBYTESZ
         byteoffset2 = 2*(plan%nx/2+1)*plan%ny*(kssta(i)-ksta)*GFLOATBYTESZ
         iret = cudaMemcpyAsyncOffHost2Dev(  plan%cu_cd_, & ! Dev
                                             byteoffset1, & ! OFFSET Dev
                                             plan%pcarr_, & ! Host
                                             byteoffset2, & ! OFFSET Host
                         plan%str_szcd_(i), pstream_(i) )
         cudaErrChk()
      END DO
      CALL GTStop(hmem); 

      CALL GTStart(hfft)
      DO i = 1,nstreams
         byteoffset1 = plan%nx*plan%ny*(kssta(i)-ksta)        *GFLOATBYTESZ
         byteoffset2 = 2*(plan%nx/2+1)*plan%ny*(kssta(i)-ksta)*GFLOATBYTESZ
         iret = GCUFFTEXECOFFC2R(plan%icuplanr_(i), plan%cu_cd_, & ! Dev
                                                    byteoffset1, & ! OFFSET
                                                    plan%cu_rd_, & ! Dev
                                                    byteoffset2)   ! OFFSET
         cudaErrChk()
      END DO
      CALL GTStop(hfft); 

      CALL GTStart(hmem)
      DO i = 1,nstreams
         byteoffset1 = plan%nx*plan%ny*(kssta(i)-ksta)*GFLOATBYTESZ
         byteoffset2 = plan%nx*plan%ny*(kssta(i)-ksta)*GFLOATBYTESZ
         iret = cudaMemCpyAsyncOffDev2Host(  plan%prarr_, & ! Host
                                             byteoffset1, & ! OFFSET Host
                                             plan%cu_rd_, & ! Dev
                                             byteoffset2, & ! OFFSET Dev
                         plan%str_szrd_(i), pstream_(i) )
         cudaErrChk()
      END DO
      DO i = 1,nstreams
         iret = cudaStreamSynchronize(pstream_(i))
      END DO
      CALL GTStop(hmem); 


      CALL GTStart(hass)
!$omp parallel do if (kend-ksta.ge.nth) private (i,j)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,plan%ny
            DO i = 1,plan%nx
              out(i,j,k) = plan%rarr(i,j,k-ksta+1)
            END DO
         END DO
      END DO
      CALL GTStop(hass)

      CALL GTStop(htot); 

      ! Update local accumulated timers:
      comtime = GTGetTime(hcom)
      ffttime = GTGetTime(hfft)
      memtime = GTGetTime(hmem)
      tratime = GTGetTime(htra)
      asstime = GTGetTime(hass)
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

      INTEGER :: ilen,jlen,klen
      INTEGER :: ilb,isize,idist
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
