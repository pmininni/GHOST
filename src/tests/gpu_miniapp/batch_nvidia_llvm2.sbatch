#!/bin/bash
#SBATCH -J gpu_miniapp_nv2
#SBATCH -p cuda
#SBATCH -n 2
#SBATCH -c 1
#SBATCH --gres=gpu:1
#SBATCH -t 00:20:00
#SBATCH -o /home/mininni/new_paradigm/GHOST/src/tests/gpu_miniapp/nv-llvm2.%j.out
cd /home/mininni/new_paradigm/GHOST/src/tests/gpu_miniapp
module purge
module load gnu15 cuda/12.4 hwloc/2.12.0
LLVM=$HOME/sw/llvm-23.1.0
OMPI=$HOME/sw/openmpi-5.0.7-cuda-llvmflang
CUDA=/opt/ohpc/pub/utils/cuda/12.4
export PATH=$OMPI/bin:$PATH
export LD_LIBRARY_PATH=$OMPI/lib:$LLVM/lib:$LLVM/lib/x86_64-unknown-linux-gnu:$CUDA/lib64:/opt/ohpc/pub/libs/gnu15/fftw/3.3.11/lib:$LD_LIBRARY_PATH
echo "===== launch geometry, spelling 1, 1 rank, N=128"
LIBOMPTARGET_INFO=-1 srun --mpi=pmix -n 1 ./build-nv-llvm/gpu_miniapp 128 2 2>&1 | grep "Launching kernel" | sort | uniq -c | sort -rn | head -12
echo "===== spelling 2 (teams loop), 2 ranks, N=128"
srun --mpi=pmix -n 2 ./build-nv-llvm-loop/gpu_miniapp 128 10
echo "===== spelling 1, 1 rank, N=128, more threads: OMP_NUM_TEAMS=80 OMP_TEAMS_THREAD_LIMIT=256"
OMP_NUM_TEAMS=80 OMP_TEAMS_THREAD_LIMIT=256 srun --mpi=pmix -n 1 ./build-nv-llvm/gpu_miniapp 128 10 | grep "device, packed\|3 kernels\|host (FFTW\|ALL TESTS"
echo "NV RUN DONE"
