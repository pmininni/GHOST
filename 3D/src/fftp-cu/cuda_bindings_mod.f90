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

 INTERFACE

!!!!!!!!!!!!!!!!!!!! Management utilities !!!!!!!!!!!!!!!!!!!!!!!!

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

!*****************************************************************
!*****************************************************************
! cudaGetDeviceProperties
!    prop : struct cudaDeviceProp pointer
!    idev : integer device id
!*****************************************************************
  INTEGER(C_INT) function cudaGetDeviceProperties(prop,idev) bind(C,name="cudaGetDeviceProperties")
    USE iso_c_binding
    USE cutypes
    IMPLICIT NONE
    INTEGER(C_INT),value :: idev
    TYPE(cudaDeviceProp) :: prop
  END FUNCTION cudaGetDeviceProperties


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


 END INTERFACE 

END MODULE cuda_bindings
