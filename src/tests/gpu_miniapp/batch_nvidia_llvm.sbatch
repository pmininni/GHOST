#!/bin/bash
#SBATCH -J gpu_miniapp_nv
#SBATCH -p cuda
#SBATCH -n 2
#SBATCH -c 1
#SBATCH --gres=gpu:1
#SBATCH -t 00:20:00
#SBATCH -o /home/mininni/new_paradigm/GHOST/src/tests/gpu_miniapp/nv-llvm.%j.out
# Runs the LLVM flang NVIDIA build of the mini-app on a Tesla P4 as a
# batch job, for when the GPU nodes have no free core for srun.
cd /home/mininni/new_paradigm/GHOST/src/tests/gpu_miniapp
module purge
module load gnu15 cuda/12.4 hwloc/2.12.0
LLVM=$HOME/sw/llvm-23.1.0
OMPI=$HOME/sw/openmpi-5.0.7-cuda-llvmflang
CUDA=/opt/ohpc/pub/utils/cuda/12.4
export PATH=$OMPI/bin:$PATH
export LD_LIBRARY_PATH=$OMPI/lib:$LLVM/lib:$LLVM/lib/x86_64-unknown-linux-gnu:$CUDA/lib64:/opt/ohpc/pub/libs/gnu15/fftw/3.3.11/lib:$LD_LIBRARY_PATH
hostname; nvidia-smi --query-gpu=name,compute_cap --format=csv,noheader
echo "===== mini-app, 1 rank, N=64"
srun --mpi=pmix -n 1 ./build-nv-llvm/gpu_miniapp 64 10
echo "===== mini-app, 2 ranks, N=64"
srun --mpi=pmix -n 2 ./build-nv-llvm/gpu_miniapp 64 10
echo "===== mini-app, 2 ranks, N=128"
srun --mpi=pmix -n 2 ./build-nv-llvm/gpu_miniapp 128 10
echo "NV RUN DONE"
