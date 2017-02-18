#include "cuda_utils.hcu"
/***************************************************************
! cuda_utils.cu
!
! CUDA interface routines. Generally, contained here are those
! methods that require at least some CUDA. Those that are
! mainly interfaces to the CUDA RTL are contained in 'cuda_bindings.f90'
!
! 2011 Duane Rosenberg & Pablo D. Mininni
!      National Center for Atmospheric Research
!      e-mail: mininni@df.uba.ar 
!
***************************************************************/
#include <cufft.h>

extern "C" {

/* Memcpy methods: */
cudaError_t cudaMemcpyHost2Dev( void *devdst, const void *hostsrc, size_t count)
{
  cudaError_t iret;
  iret = cudaMemcpy( devdst, hostsrc, count, cudaMemcpyHostToDevice ) ;
  return iret;
}

cudaError_t cudaMemcpyDev2Host( void *hostdst, const void *devsrc, size_t count)
{
  cudaError_t iret;
  iret = cudaMemcpy( hostdst , devsrc, count, cudaMemcpyDeviceToHost ); 
  return iret;
}

cudaError_t cudaMemcpyAsyncHost2Dev( void *devdst, const void *hostsrc, size_t count, cudaStream_t *stream)
{
  cudaError_t iret;
  iret = cudaMemcpyAsync( devdst, hostsrc, count, cudaMemcpyHostToDevice, *stream );
  return iret;
}

cudaError_t cudaMemcpyAsyncDev2Host( void *hostdst, const void *devsrc, size_t count, cudaStream_t *stream)
{
  cudaError_t iret;
  iret = cudaMemcpyAsync( hostdst, devsrc, count, cudaMemcpyDeviceToHost, *stream ); 
  return iret;
}

cudaError_t cudaMemcpyAsyncOffHost2Dev( void *devdst, size_t byteoffdev, const void *hostsrc, size_t byteoffhost, size_t count, cudaStream_t *stream)
{
  cudaError_t iret;
  iret = cudaMemcpyAsync( (char *) devdst + byteoffdev, (char *) hostsrc + byteoffhost, count, cudaMemcpyHostToDevice, *stream );
  return iret;
}

cudaError_t cudaMemcpyAsyncOffDev2Host( void *hostdst, size_t byteoffhost, const void *devsrc, size_t byteoffdev, size_t count, cudaStream_t *stream)
{
  cudaError_t iret;
  iret = cudaMemcpyAsync( (char *) hostdst + byteoffhost, (char *) devsrc + byteoffdev, count, cudaMemcpyDeviceToHost, *stream ); 
  return iret;
}

/* Stream methods: */
cudaError_t ptr_cudaStreamCreate( cudaStream_t **stream)
{
  *stream = (cudaStream_t *) malloc(sizeof(cudaStream_t));
  return cudaStreamCreate( *stream );
}

cudaError_t f_cudaStreamSynchronize( cudaStream_t *stream)
{
  cudaError_t iret;
  iret = cudaStreamSynchronize( *stream );
  return iret;
}

cufftResult f_cufftSetStream( cufftHandle plan, cudaStream_t *stream)
{
  cufftResult iret;
  iret = cufftSetStream( plan, *stream );
  return iret;
}

/* Interfaces for cuFFT with offsets: */

cufftResult cufftExecOffC2R( cufftHandle plan, void *datain, size_t byteoffin, void *dataout, size_t byteoffout)
{
  cufftResult iret;
  char* ptrin  = (char *) datain  + byteoffin;
  char* ptrout = (char *) dataout + byteoffout;
  iret = cufftExecC2R( plan, (cufftComplex *) ptrin, (cufftReal *) ptrout );
  return iret;
}

cufftResult cufftExecOffR2C( cufftHandle plan, void *datain, size_t byteoffin, void *dataout, size_t byteoffout)
{
  cufftResult iret;
  char* ptrin  = (char *) datain  + byteoffin;
  char* ptrout = (char *) dataout + byteoffout;
  iret = cufftExecR2C( plan, (cufftReal *) ptrin, (cufftComplex *) ptrout );
  return iret;
}

cufftResult cufftExecOffC2C( cufftHandle plan, void *datain, size_t byteoffin, void *dataout, size_t byteoffout, int dir)
{
  cufftResult iret;
  char* ptrin  = (char *) datain  + byteoffin;
  char* ptrout = (char *) dataout + byteoffout;
  iret = cufftExecC2C( plan, (cufftComplex *) ptrin, (cufftComplex *) ptrout, dir );
  return iret;
}

cufftResult cufftExecOffZ2D( cufftHandle plan, void *datain, size_t byteoffin, void *dataout, size_t byteoffout)
{
  cufftResult iret;
  char* ptrin  = (char *) datain  + byteoffin;
  char* ptrout = (char *) dataout + byteoffout;
  iret = cufftExecZ2D( plan, (cufftComplex *) ptrin, (cufftReal *) ptrout );
  return iret;
}

cufftResult cufftExecOffD2Z( cufftHandle plan, void *datain, size_t byteoffin, void *dataout, size_t byteoffout)
{
  cufftResult iret;
  char* ptrin  = (char *) datain  + byteoffin;
  char* ptrout = (char *) dataout + byteoffout;
  iret = cufftExecD2Z( plan, (cufftReal *) ptrin, (cufftComplex *) ptrout );
  return iret;
}

cufftResult cufftExecOffZ2Z( cufftHandle plan, void *datain, size_t byteoffin, void *dataout, size_t byteoffout, int dir)
{
  cufftResult iret;
  char* ptrin  = (char *) datain  + byteoffin;
  char* ptrout = (char *) dataout + byteoffout;
  iret = cufftExecZ2Z( plan, (cufftComplex *) ptrin, (cufftComplex *) ptrout, dir );
  return iret;
}

} /* end, extern "C" interface */


