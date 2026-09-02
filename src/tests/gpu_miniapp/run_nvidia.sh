#!/bin/bash
# Runs the NVIDIA build through SLURM. Usage: run_nvidia.sh [N] [niter]
# Environment: NP (ranks, default 2), NG (GPUs, default 1), NODE (g2)
cd "$(dirname "$0")"
module purge
module load nvhpc-stack/24.5
module load nvhpc-hpcx-cuda12/24.5
module load fftw/3.3.11 2>/dev/null || true
export LD_LIBRARY_PATH=/opt/ohpc/pub/libs/gnu15/fftw/3.3.11/lib:$LD_LIBRARY_PATH
NP=${NP:-2}; NG=${NG:-1}; NODE=${NODE:-g2}; BUILD=${BUILD:-build-nv}
srun --mpi=pmix -p cuda -w $NODE -n $NP -c ${CPT:-1} --gres=gpu:$NG -t 00:10:00 \
     --export=ALL ./$BUILD/gpu_miniapp "$@"
