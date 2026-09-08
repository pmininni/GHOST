#!/bin/bash
# Example that builds GHOST with OpenMP offload for AMD GPUs
# (gfx90a) using ROCm amdflang. A GPU-aware MPI library is needed.
# Assumes that MPI and ROCm compilers are in the path, and that
# mpifort points to the amdflang compiler.
# Usage: build_amd.sh [SINGLE|DOUBLE]
# P_GPU implies P_HYBRID and selects the fftp-gpu backend (hipFFT).
set -e
PREC=${1:-SINGLE}
SRC=$(cd $(dirname $0)/.. && pwd)
BLD=$SRC/build/amd
FFTDIR=/foo/fftw-3.3.10
export ROCM_PATH=/opt/rocm
mkdir -p $BLD && cd $BLD
CC=amdclang FC=mpifort cmake $SRC -DFFTW3_ROOT=$FFTDIR \
  -DP_GPU=ON -DGPU_VENDOR=AMD -DGPU_ARCH=gfx90a -DROCM_PATH=$ROCM_PATH \
  -DPRECISION=$PREC -DBIN_DIR=$SRC/bin/amd -DLIB_DIR=$SRC/lib/amd
make -j 16
