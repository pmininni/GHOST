#include <cufft.h>

/*
  cuwrappers.cu

  CUDA wrappers. This includes those structures that allow Fortran
  to access C datatypes that might vary with different version of CUDA.
  There may be some overlap between these methods, and those in 
  cuda_utils.cu.

  2013 D. Rosenberg & P. Mininni
       ORNL NCCS
       email: mininni@df.uba.edu
*/


struct cudaDevicePropG /* must agree with cudaDevicePropG in cutypes.mod */
{
  int    canMapHostMemory;
  int    clockRate;
  int    computeMode;
  int    deviceOverlap;
  int    integrated;
  int    kernelExecTimeoutEnabled;
  int    major;
  int    maxGridSize[3];
  int    maxThreadsDim[3];
  int    maxThreadsPerBlock;
  int    memoryBusWidth;
  int    memoryClockRate;
  size_t memPitch;
  int    minor;
  int    multiProcessorCount;
  char   name[256];
  int    regsPerBlock;
  size_t sharedMemPerBlock;
  size_t textureAlignment;
  size_t totalConstMem;
  size_t totalGlobalMem;
  int    warpSize;
};

extern "C" {


/*
   void w_cudagetdeviceproperties_(cudaDevicePropG *devprop, int *idev)

   Retrieves standard set of device properties

   Since cudaDeviceProp can add or delete types, we systematize them
   in the structure cudaDevicePropG, which is known both the Fortran
   and to the C:
*/
void w_cudagetdeviceproperties_(cudaDevicePropG *devprop, int *idev)
{
  cudaDeviceProp prop;

  cudaGetDeviceProperties(&prop,*idev);

  devprop->canMapHostMemory           = prop.canMapHostMemory;
  devprop->clockRate                  = prop.clockRate;
  devprop->computeMode                = prop.computeMode;
  devprop->deviceOverlap              = prop.deviceOverlap;
  devprop->integrated                 = prop.integrated;
  devprop->kernelExecTimeoutEnabled   = prop.kernelExecTimeoutEnabled;
  devprop->major                      = prop.major;
  memcpy(devprop->maxGridSize   ,prop.maxGridSize  , 3  *sizeof(int));
  memcpy(devprop->maxThreadsDim ,prop.maxThreadsDim, 3  *sizeof(int));
  devprop->maxThreadsPerBlock         = prop.maxThreadsPerBlock;
  devprop->memoryBusWidth             = prop.memoryBusWidth;
  devprop->memoryClockRate            = prop.memoryClockRate;
  devprop->memPitch                   = prop.memPitch;
  devprop->minor                      = prop.minor;
  devprop->multiProcessorCount        = prop.multiProcessorCount;
  memcpy(devprop->name          , prop.name        , 256*sizeof(char));
  devprop->regsPerBlock               = prop.regsPerBlock;
  devprop->sharedMemPerBlock          = prop.sharedMemPerBlock;
  devprop->textureAlignment           = prop.textureAlignment;
  devprop->totalConstMem              = prop.totalConstMem;
  devprop->totalGlobalMem             = prop.totalGlobalMem;
  devprop->warpSize                   = prop.warpSize;

} /* end of method w_cudagetdeviceproperties_ */


} /* extern C */

