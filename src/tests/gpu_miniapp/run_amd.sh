#!/bin/bash
# Runs the AMD build through SLURM. Usage: run_amd.sh [N] [niter]
# Environment: NP (ranks, default 2), NG (GPUs, default 1), NODE (a1),
# CPT (cores per task, 1). Keep to one node: inter-node is TCP for now.
cd "$(dirname "$0")"
module purge
module load rocm/7.2
module load openmpi5-rocm/5.0.7
export ROCM_PATH=${ROCM_PATH:-/opt/ohpc/pub/utils/rocm/7.2/rocm-7.2.3}
# The llvm/22.0 module puts the debug OpenMP runtimes (lib-debug) first
# on the library path; with them the host fallback of target regions
# segfaults. Drop that directory.
export LD_LIBRARY_PATH=$(echo $LD_LIBRARY_PATH | tr ':' '\n' | grep -v 'lib-debug' | paste -sd:)
export LD_LIBRARY_PATH=/opt/ohpc/pub/libs/gnu15/fftw/3.3.11/lib:$LD_LIBRARY_PATH
NP=${NP:-2}; NG=${NG:-1}; NODE=${NODE:-a1}; BUILD=${BUILD:-build-amd}
# GPU-aware MPI through UCX with the ROCm transports; the default ob1
# path moves device buffers about 10x slower.
export OMPI_MCA_pml=${OMPI_MCA_pml:-ucx}
export UCX_TLS=${UCX_TLS:-rocm_copy,rocm_ipc,sm,self,rc,tcp}
srun --mpi=pmix -p rocm -w $NODE -n $NP -c ${CPT:-1} --gres=gpu:$NG -t 00:10:00 \
     --export=ALL ./$BUILD/gpu_miniapp "$@"
