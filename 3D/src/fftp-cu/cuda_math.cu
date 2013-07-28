#include "cuda_math.hcu"

/***************************************************************
! cuda_math.cu
!
! CUDA math routines. Taken largely from 
!  'Optimizing Matrix Transpose in CUDA'
!  G. Ruetsch & P. Micikevicius
!  NVIDA
!
! 2013 Duane Rosenberg & Pablo D. Mininni
!      National Center for Computational Sciences
!      email: mininni@df.uba.ar 
!
***************************************************************/

__global__ void CUTranspose3(REALSZ *outdata, REALSZ *indata, int nx, int ny, int nz)
{
  __shared__ REALSZ tile[TILE_DIM][TILE_DIM][TILE_DIM+1];

  // read the matrix tile into shared memory
  unsigned int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  unsigned int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;
  unsigned int zIndex = blockIdx.z * TILE_DIM + threadIdx.z;
  if( (xIndex < nx ) && (yIndex < ny) && (zIndex < nz) ) {
    unsigned int index_in = zIndex*nx*ny + yIndex*nx + xIndex;
    tile[threadIdx.z][threadIdx.y][threadIdx.x] = indata[index_in];
  }

  __syncthreads();

  // write the transposed matrix tile to global memory
  xIndex = blockIdx.z * TILE_DIM + threadIdx.x;
  yIndex = blockIdx.y * TILE_DIM + threadIdx.y;
  zIndex = blockIdx.x * TILE_DIM + threadIdx.z;
  if( (xIndex < nz ) && (yIndex < ny) && (zIndex < nx) ) {
    unsigned int index_out = zIndex*nz*ny + yIndex*nz + xIndex;
    outdata[index_out] = tile[threadIdx.x][threadIdx.y][threadIdx.z];
  } 

} /* end of method CuTranspose3 */



__global__ void CUTranspose3C(COMPLEXSZ *outdata, COMPLEXSZ *indata, int nx, int ny, int nz)
{
  __shared__ COMPLEXSZ tile[TILE_DIM][TILE_DIM][TILE_DIM+1];

  // read the matrix tile into shared memory
  unsigned int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  unsigned int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;
  unsigned int zIndex = blockIdx.z * TILE_DIM + threadIdx.z;
  if( (xIndex < nx ) && (yIndex < ny) && (zIndex < nz) ) {
    unsigned int index_in = zIndex*nx*ny + yIndex*nx + xIndex;
    tile[threadIdx.z][threadIdx.y][threadIdx.x] = indata[index_in];
  }

  __syncthreads();

  // write the transposed matrix tile to global memory
  xIndex = blockIdx.z * TILE_DIM + threadIdx.x;
  yIndex = blockIdx.y * TILE_DIM + threadIdx.y;
  zIndex = blockIdx.x * TILE_DIM + threadIdx.z;
  if( (xIndex < nz ) && (yIndex < ny) && (zIndex < nx) ) {
    unsigned int index_out = zIndex*nz*ny + yIndex*nz + xIndex;
    outdata[index_out] = tile[threadIdx.x][threadIdx.y][threadIdx.z];
  } 

} /* end of method CuTranspose3C */


__global__ void CUTranspose(REALSZ *outdata, REALSZ *indata, int width, int height)
{
  __shared__ REALSZ tile[TILE_DIM][TILE_DIM+1];

  // read the matrix tile into shared memory
  unsigned int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  unsigned int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;
  if( (xIndex < width) && (yIndex < height) ) {
    unsigned int index_in = yIndex * width + xIndex;
    tile[threadIdx.y][threadIdx.x] = indata[index_in];
  }

  __syncthreads();

  // write the transposed matrix tile to global memory
  xIndex = blockIdx.y * TILE_DIM + threadIdx.x;
  yIndex = blockIdx.x * TILE_DIM + threadIdx.y;
  if ( (xIndex < height) && (yIndex < width) ) {
    unsigned int index_out = yIndex * height + xIndex;
    outdata[index_out] = tile[threadIdx.x][threadIdx.y];
  } 

} /* end of method CuTranspose */


