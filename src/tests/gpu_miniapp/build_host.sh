#!/bin/bash
# Host OpenMP build with gfortran + OpenMPI (validates the harness).
set -e
cd "$(dirname "$0")"
module purge
module load gnu15 openmpi5 fftw/3.3.11 cmake/4.1.2
FC=mpifort cmake -S . -B build-host -DGPU_VENDOR=NONE \
   -DFFTW3_ROOT=/opt/ohpc/pub/libs/gnu15/fftw/3.3.11 "$@"
cmake --build build-host -j 4
