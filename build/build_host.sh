#!/bin/bash
# Builds GHOST for the host with gfortran: hybrid MPI+OpenMP (default)
# or pure MPI.
# Usage: build_host.sh [hybrid|mpi] [SINGLE|DOUBLE]
# The binary goes to GHOST/bin/GHOST (BIN_DIR); the build tree to
# build/host-<mode>.
set -e
MODE=${1:-hybrid}; PREC=${2:-SINGLE}
SRC=$(cd $(dirname $0)/.. && pwd)
BLD=$SRC/build/host-$MODE
FFTDIR=/foo/fftw-3.3.10
MPIDIR=/foo/openmpi-5.0.10
[ $MODE = hybrid ] && HYB=ON || HYB=OFF
mkdir -p $BLD && cd $BLD
CC=gcc FC=gfortran cmake $SRC \
  -DFFTW3_ROOT=$FFTDIR -DCMAKE_PREFIX_PATH=$MPIDIR \
  -DP_HYBRID=$HYB -DP_GPU=OFF -DPRECISION=$PREC \
  -DBIN_DIR=$SRC/bin -DLIB_DIR=$BLD/lib
make -j 8
