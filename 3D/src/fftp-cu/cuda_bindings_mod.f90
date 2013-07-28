!=================================================================
! cuda_bindings
! Provides Fortran bindings for a variety of CUDA functions.
!
! 2011 Duane Rosenberg, Pablo D. Mininni
!      Copyright 2011-2012
!      National Center for Atmospheric Research
!      e-mail: mininni@df.uba.ar 
!
!=================================================================

MODULE cuda_bindings
 USE iso_c_binding
 IMPLICIT  NONE

 INTERFACE

!!!!!!!!!!!!!!!!!!!! Management utilities !!!!!!!!!!!!!!!!!!!!!!!!

!*****************************************************************
!*****************************************************************
! setaffinity_for_nvidia
!    my_rank: integer MPI task rank
!    ppn    : integer MPI tasks per node (MPI_PPN env. var. overrides this)
!    my_gpu : integer GPU id. IF < 0, setaffinity will choose based on
!             rank%ppn 
!*****************************************************************
  INTEGER(C_INT) function setaffinity_for_nvidia(my_rank,ppn,my_gpu) bind(C,name="setaffinity_for_nvidia")
    USE iso_c_binding
    IMPLICIT NONE
    INTEGER(C_INT)       :: my_rank
    INTEGER(C_INT)       :: ppn
    INTEGER(C_INT)       :: my_gpu
!write(*,*)'setaffinity_for_nvidia (bind): my_rank=',my_rank,' ppn=',ppn, ' my_gpu=',my_gpu
  END FUNCTION setaffinity_for_nvidia


!*****************************************************************
!*****************************************************************
! cudaGetDeviceCount
!    count : integer
!*****************************************************************
  INTEGER(C_INT) function cudaGetDeviceCount(ndev) bind(C,name="cudaGetDeviceCount")
    USE iso_c_binding
    IMPLICIT NONE
    INTEGER(C_INT)       :: ndev
  END FUNCTION cudaGetDeviceCount

!*****************************************************************
!*****************************************************************
! cudaSetDevice
!    idev : integer device id
!*****************************************************************
  INTEGER(C_INT) function cudaSetDevice(idev) bind(C,name="cudaSetDevice")
    USE iso_c_binding
    IMPLICIT NONE
    INTEGER(C_INT),value :: idev
  END FUNCTION cudaSetDevice

!*****************************************************************
!*****************************************************************
! cudaGetDevice
!    idev : integer device id
!*****************************************************************
  INTEGER(C_INT) function cudaGetDevice(idev) bind(C,name="cudaGetDevice")
    USE iso_c_binding
    IMPLICIT NONE
    INTEGER(C_INT)       :: idev
  END FUNCTION cudaGetDevice

!*****************************************************************
!*****************************************************************
! cudaSetDeviceFlags
!    flag : C_INT type
!*****************************************************************
  INTEGER(C_INT) function cudaSetDeviceFlags(flag)  bind(C,name="cudaSetDeviceFlags")
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    INTEGER(C_INT),value :: flag
  END FUNCTION cudaSetDeviceFlags

!!!*****************************************************************
!!!*****************************************************************
!   cudaGetDeviceProperties
!     gets 'standardized' set of device properties using wrapper
!     method.
!      prop : struct cudaDeviceProp pointer
!      idev : integer device id
!  *****************************************************************
    SUBROUTINE cudaGetDeviceProperties(prop,idev) bind(C,name="w_cudagetdeviceproperties_")
      USE iso_c_binding
      USE cutypes
      IMPLICIT NONE
      INTEGER(C_INT)        :: idev
      TYPE(cudaDevicePropG) :: prop
    END SUBROUTINE cudaGetDeviceProperties


!!!!!!!!!!!!!!!!!!!! Mem/compute utilities !!!!!!!!!!!!!!!!!!!!!!!

!*****************************************************************
!*****************************************************************
! cudaMemcpyHost2Dev
!    buffer: void* devdst
!    buffer: void* hstsrc
!    count : integer
!*****************************************************************
  INTEGER(C_INT) function cudaMemcpyHost2Dev(devdst, hostsrc, count)  bind(C,name="cudaMemcpyHost2Dev")
    USE iso_c_binding
    IMPLICIT NONE
    TYPE(C_PTR),value       :: devdst, hostsrc
    INTEGER(C_SIZE_T),value :: count
  END FUNCTION cudaMemcpyHost2Dev

!*****************************************************************
!*****************************************************************
!
! cudaMemcpyDev2Host
!    buffer: void* devdst
!    buffer: void* hstsrc
!    count : integer
!*****************************************************************
  INTEGER(C_INT) function cudaMemcpyDev2Host(hostdst, devsrc, count)  bind(C,name="cudaMemcpyDev2Host")
    USE iso_c_binding
    IMPLICIT NONE
    TYPE(C_PTR),value      :: devsrc, hostdst
    INTEGER(C_SIZE_T),value:: count
  END FUNCTION cudaMemcpyDev2Host
 
!*****************************************************************
!*****************************************************************
! cudaMalloc
!    buffer: void** pointer
!    isize : integer byte size of buffer 
!*****************************************************************
  INTEGER(C_INT) function cudaMalloc(buffer, isize)  bind(C,name="cudaMalloc")
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    TYPE(C_PTR)              :: buffer
    INTEGER(C_SIZE_T), value :: isize
  END FUNCTION cudaMalloc
 
!*****************************************************************
!*****************************************************************
! cudaHostAlloc
!    buffer: void** pointer
!    isize : integer byte size of buffer
!    flag  : one of cudaHostAllocMapped (recommended), cudaHostAllocWriteCombined,
!            cudaHostAllocPortable, cudaHostAllocDefault. If the last
!            is used, this function seems to behave the same as cudaMallocHost.
!*****************************************************************
  INTEGER(C_INT) function cudaHostAlloc(buffer, isize, flags)  bind(C,name="cudaHostAlloc")
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    TYPE(C_PTR)            :: buffer
    INTEGER(C_SIZE_T),value:: isize
    INTEGER(C_INT),value   :: flags
  END FUNCTION cudaHostAlloc

!*****************************************************************
!*****************************************************************
! cudaMallocHost
!    buffer: void** pointer
!    isize : integer byte size of buffer
!*****************************************************************
  INTEGER(C_INT) function cudaMallocHost(buffer, isize)  bind(C,name="cudaMallocHost")
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    TYPE(C_PTR)               :: buffer
    INTEGER(C_SIZE_T), value  :: isize
  END FUNCTION cudaMallocHost

 
!*****************************************************************
!*****************************************************************
! cudaFreeHost
!    buffer: void* pointer
!*****************************************************************
  INTEGER(C_INT) function cudaFreeHost(buffer)  bind(C,name="cudaFreeHost")
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    TYPE(C_PTR),value :: buffer
  END FUNCTION cudaFreeHost

!*****************************************************************
!*****************************************************************
! cudaFree
!    buffer: void* pointer
!*****************************************************************
  INTEGER(C_INT) function cudaFree(buffer)  bind(C,name="cudaFree")
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    TYPE (C_PTR),value :: buffer
  END FUNCTION cudaFree

!*****************************************************************
!*****************************************************************
! cudaHostGetDevicePointer
!    pdevice: void** device pointer
!    phost  : void*  host pointer mapping
!    flags  : 0 (for now)
!*****************************************************************
  INTEGER(C_INT) function cudaHostGetDevicePointer(pdevice, phost, flags)  bind(C,name="cudaHostGetDevicePointer")
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    TYPE(C_PTR)          :: pdevice, phost
    INTEGER(C_INT), value:: flags
  END FUNCTION cudaHostGetDevicePointer

!*****************************************************************
!*****************************************************************
! cufftPlan1d :: create CUDA FFT 1d plan
!     pplan    : cufftHandle pointer (out)
!     nx       : transform size
!     cutype   : transform data type: CUFFT_R2C (0x2a), CUFFT_C2R(0x2c), CUFFT_C2C (0x29)
!     batch    : number of transforms of size nx
!
!*****************************************************************
 INTEGER(C_INT) function cufftPlan1d(pplan, nx, cutype, batch) bind(C,name="cufftPlan1d")
   USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE
  INTEGER(C_INT),value :: nx
  INTEGER(C_INT),value :: cutype
  INTEGER(C_INT),value :: batch
  INTEGER(C_INT)       :: pplan
 END FUNCTION cufftPlan1d

!*****************************************************************
!*****************************************************************
! cufftPlan2d :: create CUDA FFT 2d plan
!     pplan    : cufftHandle pointer (out)
!     nx       : transform size
!     ny       : transform size
!     cutype   : transform data type: CUFFT_R2C (0x2a), CUFFT_C2R(0x2c), CUFFT_C2C (0x29)
!
!*****************************************************************
 INTEGER(C_INT) function cufftPlan2d(pplan, nx, ny, cutype) bind(C,name="cufftPlan2d")
   USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE
  INTEGER(C_INT),value :: nx
  INTEGER(C_INT),value :: ny
  INTEGER(C_INT),value :: cutype
  INTEGER(C_INT)       :: pplan
 END FUNCTION cufftPlan2d

!*****************************************************************
!*****************************************************************
! cufftPlanMany :: create CUDA FFT 2d plan
!     pplan      : cufftHandle pointer (out)
!     rank       : integer rank
!     nx         : transform size
!     ny         : transform size
!     cutype     : transform data type: CUFFT_R2C (0x2a), CUFFT_C2R(0x2c), CUFFT_C2C (0x29)
!
!*****************************************************************
  INTEGER(C_INT) function cufftPlanMany(pplan,rank,pn,&
 pinembed,istride,idist,ponembed,ostride,odist,cutype,batch)  bind(C,name="cufftPlanMany")
!  USE, INTRINSIC :: iso_c_binding
   USE  :: iso_c_binding
  IMPLICIT NONE
  INTEGER(C_INT)        :: pplan
  INTEGER(C_INT),value  :: rank
  INTEGER(C_INT)        :: pn(*)
  INTEGER(C_INT)        :: pinembed(*)
  INTEGER(C_INT),value  :: istride
  INTEGER(C_INT),value  :: idist
  INTEGER(C_INT)        :: ponembed(*)
  INTEGER(C_INT),value  :: ostride
  INTEGER(C_INT),value  :: odist
  INTEGER(C_INT),value  :: cutype
  INTEGER(C_INT),value  :: batch
  END FUNCTION cufftPlanMany

!*****************************************************************
!*****************************************************************
! cufftDestroy:: destroy CUDA FFT plan
!     pplan     : cufftHandle pointer (in)
!
!*****************************************************************
 INTEGER(C_INT) function cufftDestroy(pplan) bind(C,name="cufftDestroy")
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE
  INTEGER(C_INT),value :: pplan
 END FUNCTION cufftDestroy

!*****************************************************************
!*****************************************************************
! cufftExecR2C:: execute CUDA FFT plan, R-->C
!    pplan     : cufftHandle pointer (out)
!    datain    : input data (cufftReal)
!    dataout   : output data (cufftComplex)
!
!*****************************************************************
 INTEGER(C_INT) function cufftExecR2C(pplan, datain, dataout) bind(C,name="cufftExecR2C")
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE
  TYPE(C_PTR),value   :: datain
  TYPE(C_PTR),value   :: dataout
  INTEGER(C_INT),value:: pplan
 END FUNCTION cufftExecR2C

!*****************************************************************
!*****************************************************************
! cufftExecC2R:: execute CUDA FFT plan, C-->R
!    pplan     : cufftHandle pointer (out)
!    datain    : input data (cufftReal)
!    dataout   : output data (cufftComplex)
!
!*****************************************************************
 INTEGER(C_INT) function cufftExecC2R(pplan, datain, dataout) bind(C,name="cufftExecC2R")
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE
  TYPE(C_PTR),value    :: datain
  TYPE(C_PTR),value    :: dataout
  INTEGER(C_INT),value :: pplan
 END FUNCTION cufftExecC2R

!*****************************************************************
!*****************************************************************
! cufftExecC2C:: execute CUDA FFT plan, C-->C
!    pplan     : cufftHandle pointer (out)
!    datain    : input data (cufftReal)
!    dataout   : output data (cufftComplex)
!    idir      : transpform direction: CUFFT_FORWARD or CUFFT_INVERSE, 
!                same as FFTCU_REAL_TO_COMPLEX, and FFTCU_COMPLEX_TO_REAL, respectively
!
!*****************************************************************
 INTEGER(C_INT) function cufftExecC2C(pplan, datain, dataout, idir) bind(C,name="cufftExecC2C")
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE
  TYPE(C_PTR),value    :: datain
  TYPE(C_PTR),value    :: dataout
  INTEGER(C_INT),value :: idir
  INTEGER(C_INT),value :: pplan
 END FUNCTION cufftExecC2C

!*****************************************************************
!*****************************************************************
! cufftExecD2Z:: execute CUDA FFT plan, D-->Z, double to double complex
!    pplan     : cufftHandle pointer (out)
!    datain    : input data (cufftDouble)
!    dataout   : output data (cufftDoubleComplex)
!
!*****************************************************************
 INTEGER(C_INT) function cufftExecD2Z(pplan, datain, dataout) bind(C,name="cufftExecD2Z")
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE
  TYPE(C_PTR),value   :: datain
  TYPE(C_PTR),value   :: dataout
  INTEGER(C_INT),value:: pplan
 END FUNCTION cufftExecD2Z

!*****************************************************************
!*****************************************************************
! cufftExecZ2D:: execute CUDA FFT plan, Z-->D, double complex to double
!    pplan     : cufftHandle pointer (out)
!    datain    : input data (cufftDoubleComplex)
!    dataout   : output data (cufftDouble)
!
!*****************************************************************
 INTEGER(C_INT) function cufftExecZ2D(pplan, datain, dataout) bind(C,name="cufftExecZ2D")
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE
  TYPE(C_PTR),value    :: datain
  TYPE(C_PTR),value    :: dataout
  INTEGER(C_INT),value :: pplan
 END FUNCTION cufftExecZ2D

!*****************************************************************
!*****************************************************************
! cufftExecZ2Z:: execute CUDA FFT plan, Z-->Z, double complex to double complex
!    pplan     : cufftHandle pointer (out)
!    datain    : input data (cufftDoubleComplex)
!    dataout   : output data (cufftDoubleComplex)
!    idir      : transpform direction: CUFFT_FORWARD or CUFFT_INVERSE, 
!                same as FFTCU_REAL_TO_COMPLEX, and FFTCU_COMPLEX_TO_REAL, respectively
!
!*****************************************************************
 INTEGER(C_INT) function cufftExecZ2Z(pplan, datain, dataout, idir) bind(C,name="cufftExecZ2Z")
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE
  TYPE(C_PTR),value    :: datain
  TYPE(C_PTR),value    :: dataout
  INTEGER(C_INT),value :: idir
  INTEGER(C_INT),value :: pplan
 END FUNCTION cufftExecZ2Z


!*****************************************************************
!*****************************************************************
! cuTranspose :: take tranpose of 'datain' and store in 'dataout'. Data
!                must reside on device.
!    dataout   : output data 
!    datain    : input data 
!    width     : no columns of input matrix
!    height    : no. rows  of input matrix
!
!*****************************************************************
 SUBROUTINE cuTranspose(dataout, datain, width, height) bind(C,name="w_cudatranspose_")
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE
  TYPE(C_PTR),value    :: dataout
  TYPE(C_PTR),value    :: datain
  INTEGER(C_INT),value :: width
  INTEGER(C_INT),value :: height
 END SUBROUTINE cuTranspose


!*****************************************************************
!*****************************************************************
! cuTranspose3:: take tranpose of real 'datain' and store in 'dataout',
!                transposing x <-> z. Data must reside on device.
!    dataout   : output data 
!    datain    : input data 
!    nx,ny,nz  : dimensions of input array
!
!*****************************************************************
 SUBROUTINE cuTranspose3(dataout, datain, nx, ny, nz) bind(C,name="w_cudatranspose3_")
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE
  TYPE(C_PTR),value    :: dataout
  TYPE(C_PTR),value    :: datain
  INTEGER(C_INT),value :: nx, ny, nz
 END SUBROUTINE cuTranspose3


!*****************************************************************
!*****************************************************************
! cuTranspose3C:: take tranpose of complex 'datain' and store in 'dataout',
!                 transposing x <-> z. Data must reside on device.
!    dataout   : output data 
!    datain    : input data 
!    nx,ny,nz  : dimensions of input array
!
!*****************************************************************
 SUBROUTINE cuTranspose3C(dataout, datain, nx, ny, nz) bind(C,name="w_cudatranspose3C_")
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE
  TYPE(C_PTR),value    :: dataout
  TYPE(C_PTR),value    :: datain
  INTEGER(C_INT),value :: nx, ny, nz
 END SUBROUTINE cuTranspose3C


 END INTERFACE 


 CONTAINS


!*****************************************************************
!*****************************************************************
  SUBROUTINE GetCUFFTErr(iret, sret)
!    PURPOSE: Get error string condition for CUFFT calls.
!    iret : integer return value from CUFFT function (IN)
!    sret : string of size >= 24 to hold error condition (OUT)
!*****************************************************************
      INTEGER, INTENT(IN)                     :: iret
      CHARACTER(len=*), INTENT(OUT)           :: sret

      CHARACTER(len=24), DIMENSION(0:8)       :: serr

      serr(0) = 'CUFFT_SUCCESS'
      serr(1) = 'CUFFT_INVALID_PLAN'
      serr(2) = 'CUFFT_ALLOC_FAILED'
      serr(3) = 'CUFFT_INVALID_TYPE'
      serr(4) = 'CUFFT_INVALID_VALUE'
      serr(5) = 'CUFFT_INTERNAL_ERROR'
      serr(6) = 'CUFFT_EXEC_FAILED CUFFT'
      serr(7) = 'CUFFT_SETUP_FAILED'
      serr(8) = 'CUFFT_INVALID_SIZE'

      IF ( iret .LT. 0 .OR. iret .GT. 8 ) THEN
        sret = 'UNKNOWN_ERROR'
        RETURN
      ENDIF

      sret = serr(iret)

  END SUBROUTINE GetCUFFTErr


END MODULE cuda_bindings
