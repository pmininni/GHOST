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
cudaError_t cudaMemcpyHost2Dev( void *devdst, const void *hostsrc,  size_t count)
{
  cudaError_t iret;
  iret = cudaMemcpy( devdst, hostsrc, count,  cudaMemcpyHostToDevice ) ;
  return iret;
}

cudaError_t cudaMemcpyDev2Host( void *hostdst, const void *devsrc,  size_t count)
{
  cudaError_t iret;
  iret = cudaMemcpy( hostdst , devsrc, count,  cudaMemcpyDeviceToHost ); 
  return iret;
}

/* New interface for cufftPlanMany: */
int  MycufftPlanMany(cufftHandle *plan, int rank, int *n, int *inembed,
                       int istride, int idist, int *onembed, int ostride,
                       int odist, cufftType type, int batch) 
{
  int iret;
  
  if ( inembed != NULL && inembed[0] <= 0 ) {
    iret = (int)cufftPlanMany(plan,rank,n,NULL,istride,idist,
                                     NULL,ostride,odist,type,batch);
  } else {
    iret = (int)cufftPlanMany(plan,rank,n,inembed,istride,idist,
                                           onembed,ostride,odist,type,batch);
  }
  return iret;
}


} /* end, extern "C" interface */


