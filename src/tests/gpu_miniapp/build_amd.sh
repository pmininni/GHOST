#!/bin/bash
# AMD offload build: ROCm 7.2 amdflang with a ROCm-aware OpenMPI that
# was itself built with amdflang (the system one was built with
# gfortran, whose mpi_f08 module amdflang cannot read).
set -e
cd "$(dirname "$0")"
module purge
module load rocm/7.2 cmake/4.1.2 hwloc/2.12.0
export ROCM_PATH=${ROCM_PATH:-/opt/ohpc/pub/utils/rocm/7.2/rocm-7.2.3}
OMPI=${OMPI:-$HOME/sw/openmpi-5.0.7-rocm-amdflang}
export PATH=$OMPI/bin:$ROCM_PATH/lib/llvm/bin:$PATH
export LD_LIBRARY_PATH=$OMPI/lib:$LD_LIBRARY_PATH
BUILD=${BUILD:-build-amd}
FC=mpifort cmake -S . -B $BUILD -DGPU_VENDOR=AMD -DGPU_ARCH=${GPU_ARCH:-gfx90a} \
   -DROCM_PATH=$ROCM_PATH \
   -DFFTW3_ROOT=/opt/ohpc/pub/libs/gnu15/fftw/3.3.11 "$@"
cmake --build $BUILD -j 4
