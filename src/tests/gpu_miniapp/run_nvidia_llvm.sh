#!/bin/bash
# Runs the LLVM flang NVIDIA build through SLURM. Usage: run_nvidia_llvm.sh [N] [niter]
# Environment: NP (ranks, default 2), NG (GPUs, default 1), NODE (g2), CPT (cores/task, 1)
cd "$(dirname "$0")"
module purge
module load gnu15 cuda/12.4 hwloc/2.12.0
LLVM=${LLVM:-$HOME/sw/llvm-23.1.0}
OMPI=${OMPI:-$HOME/sw/openmpi-5.0.7-cuda-llvmflang}
CUDA=${CUDA_HOME:-/opt/ohpc/pub/utils/cuda/12.4}
export PATH=$OMPI/bin:$PATH
export LD_LIBRARY_PATH=$OMPI/lib:$LLVM/lib:$LLVM/lib/x86_64-unknown-linux-gnu:$CUDA/lib64:/opt/ohpc/pub/libs/gnu15/fftw/3.3.11/lib:$LD_LIBRARY_PATH
NP=${NP:-2}; NG=${NG:-1}; NODE=${NODE:-g2}; BUILD=${BUILD:-build-nv-llvm}
srun --mpi=pmix -p cuda -w $NODE -n $NP -c ${CPT:-1} --gres=gpu:$NG -t 00:10:00 \
     --export=ALL ./$BUILD/gpu_miniapp "$@"
