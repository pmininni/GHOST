#!/bin/bash
# Example that builds GHOST with OpenMP offload for the AMD GPUs 
# (gfx90a) with ROCm amdflang. A GPU-aware MPI library is needed.
# Usage: build_amd.sh [SINGLE|DOUBLE]
set -e
PREC=${1:-SINGLE}
SRC=$(cd $(dirname $0)/.. && pwd)
BLD=$SRC/build/amd
module purge; module load rocm/7.2
module load llvm/22.0 openmpi5-rocm/5.0.7 cmake/4.1.2 fftw-gnu/3.3.11
export ROCM_PATH=/opt/ohpc/pub/utils/rocm/7.2/rocm-7.2.3
mkdir -p $BLD && cd $BLD
CC=amdclang FC=mpifort cmake $SRC \
  -DFFTW3_ROOT=/opt/ohpc/pub/libs/gnu15/fftw/3.3.11 \
  -DP_GPU=ON -DGPU_VENDOR=AMD -DGPU_ARCH=gfx90a -DROCM_PATH=$ROCM_PATH \
  -DPRECISION=$PREC -DBIN_DIR=$SRC/bin/amd -DLIB_DIR=$SRC/lib/amd
make -j 16
