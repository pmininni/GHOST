#!/bin/bash
# Example that builds GHOST with OpenMP offload for NVIDIA GPUs.
# Assumes that MPI is in the path, that mpifort points to the right
# compiler, and that LLVM or nvfortran compilers are in the path.
# Usage: build_nvidia.sh [llvm|nvhpc] [ARCH] [SINGLE|DOUBLE]
#   llvm : LLVM flang with NVPTX support and a CUDA-aware MPI library.
#   nvhpc: nvfortran with its own MPI; needs compute capability >=7.0
# P_GPU implies P_HYBRID and selects the fftp-gpu backend (cuFFT).
set -e
TC=${1:-llvm}; PREC=${3:-SINGLE}
SRC=$(cd $(dirname $0)/.. && pwd)
BLD=$SRC/build/nvidia-$TC
mkdir -p $BLD && cd $BLD
if [ $TC = llvm ]; then
  ARCH=${2:-sm_61}
  FFTDIR=/foo/fftw-3.3.10
  CUDA_HOME=/foo/cuda-12.4
  CC=clang FC=mpifort cmake $SRC -DFFTW3_ROOT=$FFTDIR \
    -DP_GPU=ON -DGPU_VENDOR=NVIDIA -DGPU_ARCH=$ARCH -DCUDA_PATH=$CUDA_HOME \
    -DPRECISION=$PREC -DBIN_DIR=$SRC/bin/nvidia -DLIB_DIR=$SRC/lib/nvidia
else
  ARCH=${2:-cc70}
  FFTDIR=/foo/fftw-3.3.10
  # NVHPC compilers older than 24.11 have a bug in use_device_addr
  # and don't offload properly, use LLVM compilers instead if you 
  # need support for older GPUs or CUDA versions.
  NV=/foo/nvidia/hpc_sdk/Linux_x86_64/24.11
  CC=nvc FC=$NV/comm_libs/12.6/hpcx/ompi/bin/mpifort cmake $SRC \
    -DFFTW3_ROOT=$FFTDIR -DP_GPU=ON -DGPU_VENDOR=NVIDIA -DGPU_ARCH=$ARCH \
    -DCUDA_PATH=$NV/cuda/12.6 -DPRECISION=$PREC \
    -DBIN_DIR=$SRC/bin/nvidia -DLIB_DIR=$SRC/lib/nvidia
fi
make -j 16
