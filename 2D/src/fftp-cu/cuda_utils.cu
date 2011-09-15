#include "cuda_utils.hcu"

extern "C" {

/* Memcpy methods: */
cudaError_t cudaMemcpyHost2Dev( void *devdst, const void *hostsrc,  size_t count)
{
  cudaError_t iret;
  iret = cudaMemcpy( devdst, hostsrc, count,  cudaMemcpyHostToDevice ) ;
  return iret;
}

cudaError_t cudaMemcpyDev2Host( void *hostdst, const void *devsrc,  size_t count)
{
  return cudaMemcpy( hostdst , devsrc, count,  cudaMemcpyDeviceToHost ); 
}


} /* end, extern "C" interface */


