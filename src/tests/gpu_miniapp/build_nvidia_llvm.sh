#!/bin/bash
# NVIDIA offload build with upstream LLVM flang (NVPTX backend), the only
# OpenMP route to the Tesla P4 (cc 6.1): nvfortran refuses OpenMP offload
# below cc70. Uses a CUDA-aware OpenMPI built with the same flang.
set -e
cd "$(dirname "$0")"
module purge
module load gnu15 cmake/4.1.2 cuda/12.4 hwloc/2.12.0
LLVM=${LLVM:-$HOME/sw/llvm-23.1.0}
OMPI=${OMPI:-$HOME/sw/openmpi-5.0.7-cuda-llvmflang}
export CUDA_PATH=${CUDA_HOME:-/opt/ohpc/pub/utils/cuda/12.4}
export PATH=$OMPI/bin:$LLVM/bin:$PATH
export LD_LIBRARY_PATH=$OMPI/lib:$LLVM/lib:$LD_LIBRARY_PATH
BUILD=${BUILD:-build-nv-llvm}
FC=mpifort cmake -S . -B $BUILD -DGPU_VENDOR=NVIDIA -DGPU_ARCH=${GPU_ARCH:-sm_61} \
   -DCUDA_PATH=$CUDA_PATH \
   -DFFTW3_ROOT=/opt/ohpc/pub/libs/gnu15/fftw/3.3.11 "$@"
cmake --build $BUILD -j 4
