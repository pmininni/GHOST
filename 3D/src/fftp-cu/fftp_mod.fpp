!=================================================================
! MODULES for FFTP v3
! Parallel Fast Fourier Transform in 3D using CUDA
!
! 2011 Duane L. Rosenberg and Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.uba.ar 
!=================================================================

!=================================================================

  MODULE fftplans
!
! Set the variable ikind to:  4 in 32 bits machines
!                             8 in 64 bits machines
! Set the variable csize to:  8 if L1 cache is <= 64 kb
!                            16 if L1 cache is 128 kb
! The variable nstrip controls strip mining during the 
! transposition. Often set to 1.
!
      USE fprecision
      USE iso_c_binding
 
      INTEGER, PARAMETER  :: ikind = IKIND_
      INTEGER, PARAMETER  :: csize = CSIZE_
      INTEGER, PARAMETER  :: nstrip = NSTRIP_
      INTEGER, PARAMETER  :: nstreams = NSTREAMS_
      INTEGER, PARAMETER  :: FFTCU_REAL_TO_COMPLEX = -1
      INTEGER, PARAMETER  :: FFTCU_COMPLEX_TO_REAL =  1
      INTEGER, PARAMETER  :: FFTW_REAL_TO_COMPLEX = -1
      INTEGER, PARAMETER  :: FFTW_COMPLEX_TO_REAL =  1
      INTEGER, PARAMETER  :: FFTW_MEASURE =  0
      INTEGER, PARAMETER  :: FFTW_PATIENT =  0
      INTEGER, PARAMETER  :: FFTW_ESTIMATE=  0
      INTEGER, SAVE       :: streams_created = 0
      TYPE  (C_PTR)       :: pstream_(nstreams)
      INTEGER, DIMENSION (nstreams) :: issta,issnd
      INTEGER, DIMENSION (nstreams) :: kssta,kssnd
      INTEGER             :: hcom,hfft,hmem,htra,hass,htot
      DOUBLE PRECISION    :: comtime = 0.0
      DOUBLE PRECISION    :: ffttime = 0.0
      DOUBLE PRECISION    :: memtime = 0.0
      DOUBLE PRECISION    :: tratime = 0.0
      DOUBLE PRECISION    :: asstime = 0.0
      DOUBLE PRECISION    :: tottime = 0.0
      CHARACTER, TARGET   :: cfileerr(1024)
      TYPE FFTPLAN
         COMPLEX(KIND=GP), POINTER, DIMENSION (:,:,:)  :: ccarr
         COMPLEX(KIND=GP), POINTER, DIMENSION (:,:,:)  :: ccarrt
         COMPLEX(KIND=GP), POINTER, DIMENSION (:,:,:)  :: carr
         REAL   (KIND=GP), POINTER, DIMENSION (:,:,:)  :: rarr
         TYPE     (C_PTR)                              :: cu_ccd_,cu_ccd1_
         TYPE     (C_PTR)                              :: cu_cd_,cu_rd_
         TYPE     (C_PTR)                              :: pccarr_,pccarrt_,pcarr_
         TYPE     (C_PTR)                              :: prarr_
         INTEGER  (C_INT),        DIMENSION (nstreams) :: icuplanr_
         INTEGER  (C_INT),        DIMENSION (nstreams) :: icuplanc_
         INTEGER                                       :: nx,ny,nz
         INTEGER(C_SIZE_T)                             :: szccd_,szcd_,szrd_
         INTEGER(C_SIZE_T),       DIMENSION (nstreams) :: str_szccd_
         INTEGER(C_SIZE_T),       DIMENSION (nstreams) :: str_szcd_
         INTEGER(C_SIZE_T),       DIMENSION (nstreams) :: str_szrd_
         INTEGER, DIMENSION (:), POINTER               :: itype1,itype2
      END TYPE FFTPLAN
      SAVE

  END MODULE fftplans
!=================================================================

  MODULE cutypes
      USE ISO_C_BINDING
      IMPLICIT NONE
      
   
      INTEGER, PARAMETER :: CUFFT_SUCCESS       =X'0'
      INTEGER, PARAMETER :: CUFFT_INVALID_PLAN  =X'1'
      INTEGER, PARAMETER :: CUFFT_ALLOC_FAILED  =X'2'
      INTEGER, PARAMETER :: CUFFT_INVALID_TYPE  =X'3'
      INTEGER, PARAMETER :: CUFFT_INVALID_VALUE =X'4'
      INTEGER, PARAMETER :: CUFFT_INTERNAL_ERROR=X'5'
      INTEGER, PARAMETER :: CUFFT_EXEC_FAILED   =X'6'
      INTEGER, PARAMETER :: CUFFT_SETUP_FAILED  =X'7'
      INTEGER, PARAMETER :: CUFFT_INVALID_SIZE  =X'8'
      INTEGER, PARAMETER :: CUFFT_UNALIGNED_DATA=X'9'

      ! 
      INTEGER, PARAMETER :: CUFFT_R2C=X'2a'
      INTEGER, PARAMETER :: CUFFT_C2R=X'2c'
      INTEGER, PARAMETER :: CUFFT_C2C=X'29'
      INTEGER, PARAMETER :: CUFFT_D2Z=X'6a'
      INTEGER, PARAMETER :: CUFFT_Z2D=X'6c'
      INTEGER, PARAMETER :: CUFFT_Z2Z=X'69'
      INTEGER, PARAMETER :: cudaHostAllocDefault =0
      INTEGER, PARAMETER :: cudaHostAllocPortable=1
      INTEGER, PARAMETER :: cudaHostAllocMapped  =2
      INTEGER, PARAMETER :: cudaDeviceMapHost    =3
     
      ENUM, BIND(C)
        ENUMERATOR ::                            &
        cudaSuccess                       =0 ,   &
        cudaErrorMissingConfiguration     =1 ,   &
        cudaErrorMemoryAllocation         =2 ,   &
        cudaErrorInitializationError      =3 ,   &
        cudaErrorLaunchFailure            =4 ,   &
        cudaErrorPriorLaunchFailure       =5 ,   &
        cudaErrorLaunchTimeout            =6 ,   &
        cudaErrorLaunchOutOfResources     =7 ,   &
        cudaErrorInvalidDeviceFunction    =8 ,   &
        cudaErrorInvalidConfiguration     =9 ,   &
        cudaErrorInvalidDevice            =10,   &
        cudaErrorInvalidValue             =11,   &
        cudaErrorInvalidPitchValue        =12,   &
        cudaErrorInvalidSymbol            =13,   &
        cudaErrorMapBufferObjectFailed    =14,   &
        cudaErrorUnmapBufferObjectFailed  =15,   &
        cudaErrorInvalidHostPointer       =16,   &
        cudaErrorInvalidDevicePointer     =17,   &
        cudaErrorInvalidTexture           =18,   &
        cudaErrorInvalidTextureBinding    =19,   &
        cudaErrorInvalidChannelDescriptor =20,   &
        cudaErrorInvalidMemcpyDirection   =21,   &
        cudaErrorAddressOfConstant        =22,   &
        cudaErrorTextureFetchFailed       =23,   &
        cudaErrorTextureNotBound          =24,   &
        cudaErrorSynchronizationError     =25,   &
        cudaErrorInvalidFilterSetting     =26,   &
        cudaErrorInvalidNormSetting       =27,   &
        cudaErrorMixedDeviceExecution     =28,   &
        cudaErrorCudartUnloading          =29,   &
        cudaErrorUnknown                  =30,   &
        cudaErrorNotYetImplemented        =31,   &
        cudaErrorMemoryValueTooLarge      =32,   &
        cudaErrorInvalidResourceHandle    =33,   &
        cudaErrorNotReady                 =34,   &
        cudaErrorInsufficientDriver       =35,   &
        cudaErrorSetOnActiveProcess       =36,   &
        cudaErrorNoDevice                 =37,   &
        cudaErrorStartupFailure           =38,   &
        cudaErrorApiFailureBase           =39
      END ENUM

      TYPE, BIND(C) :: cudaDevicePropG
        INTEGER   (C_INT) :: canMapHostMemory
        INTEGER   (C_INT) :: managedMemory
        INTEGER   (C_INT) :: concurrentManagedAccess
        INTEGER   (C_INT) :: clockRate
        INTEGER   (C_INT) :: computeMode
        INTEGER   (C_INT) :: deviceOverlap
        INTEGER   (C_INT) :: integrated
        INTEGER   (C_INT) :: kernelExecTimeoutEnabled
        INTEGER   (C_INT) :: major
        INTEGER   (C_INT) :: maxGridSize(3)
        INTEGER   (C_INT) :: maxThreadsDim(3)
        INTEGER   (C_INT) :: maxThreadsPerBlock
        INTEGER(C_SIZE_T) :: memPitch
        INTEGER   (C_INT) :: minor
        INTEGER   (C_INT) :: multiProcessorCount
        CHARACTER(C_CHAR) :: name(256) 
        INTEGER   (C_INT) :: regsPerBlock
        INTEGER(C_SIZE_T) :: sharedMemPerBlock
        INTEGER(C_SIZE_T) :: textureAlignment
        INTEGER(C_SIZE_T) :: totalConstMem
        INTEGER(C_SIZE_T) :: totalGlobalMem
        INTEGER   (C_INT) :: warpSize
      END TYPE cudaDevicePropG

  END MODULE cutypes
!=================================================================

  MODULE threads
!
! nth: number of threads used in OpenMP loops and FFTs
!$    USE omp_lib
      INTEGER :: nth
      SAVE

  END MODULE threads

!=================================================================

  MODULE offloading
!
! Module with variables for offloading to GPUs using OpenMP. GHOST
! assumes the number of MPI jobs in each node is equal to the
! number of GPUs available in the node. The user must ensure this
! condition is fulfilled.
!$    USE threads
      INTEGER                            :: numdev
      INTEGER                            :: hostdev,targetdev
      LOGICAL                            :: offload
      SAVE
    CONTAINS
      SUBROUTINE init_offload(myrank,ndev,hdev,tdev)
      ! Binds each MPI job to only one GPU target device
        INTEGER, INTENT(IN)    :: myrank
        INTEGER, INTENT(INOUT) :: ndev,hdev,tdev
        ndev = 0
        hdev = 0
        tdev = 0
#if defined(DO_HYBRIDoffl)
!$      ndev = omp_get_num_devices()
!$      hdev = omp_get_initial_device()
!$      tdev = MODULO(myrank,ndev) + 1
#endif
      END SUBROUTINE init_offload

  END MODULE offloading
!=================================================================

  MODULE mpivars
!     INCLUDE 'mpif.h'
      INTEGER, SAVE :: ista,iend
      INTEGER, SAVE :: jsta,jend
      INTEGER, SAVE :: ksta,kend
      INTEGER, SAVE :: nprocs,myrank
      INTEGER, SAVE :: provided
      INTEGER, SAVE :: ierr

  END MODULE mpivars

!=================================================================

