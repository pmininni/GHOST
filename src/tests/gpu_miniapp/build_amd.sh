#!/bin/bash
# AMD offload build: ROCm 7.2 amdflang with the system ROCm-aware
# OpenMPI (module openmpi5-rocm, built with amdflang).
set -e
cd "$(dirname "$0")"
module purge
module load rocm/7.2 cmake/4.1.2
module load openmpi5-rocm/5.0.7
export ROCM_PATH=${ROCM_PATH:-/opt/ohpc/pub/utils/rocm/7.2/rocm-7.2.3}
BUILD=${BUILD:-build-amd}
FC=mpifort cmake -S . -B $BUILD -DGPU_VENDOR=AMD -DGPU_ARCH=${GPU_ARCH:-gfx90a} \
   -DROCM_PATH=$ROCM_PATH \
   -DFFTW3_ROOT=/opt/ohpc/pub/libs/gnu15/fftw/3.3.11 "$@"
cmake --build $BUILD -j 4
