#!/bin/bash
# NVIDIA offload build: nvfortran 24.5 with the CUDA-aware HPC-X MPI.
# Extra cmake options can be passed as arguments, e.g.
#   ./build_nvidia.sh -DGPU_LOOP_SPELLING=2 -DPRECISION=DOUBLE
set -e
cd "$(dirname "$0")"
module purge
module load nvhpc-stack/24.5
module load nvhpc-hpcx-cuda12/24.5
module load cmake/4.1.2
BUILD=${BUILD:-build-nv}
FC=mpifort cmake -S . -B $BUILD -DGPU_VENDOR=NVIDIA -DGPU_ARCH=${GPU_ARCH:-cc60} \
   -DFFTW3_ROOT=/opt/ohpc/pub/libs/gnu15/fftw/3.3.11 "$@"
cmake --build $BUILD -j 4
