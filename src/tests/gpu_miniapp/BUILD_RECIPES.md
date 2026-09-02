# Build recipes for the GPU toolchains used by the mini-app

Exact commands used on sakura, in the order they were run. Paths under
`$HOME/sw` are the private installs; replace the prefixes for a
system-wide install. Sources are in `$HOME/sw/src`.

## 1. OpenMPI 5.0.7 for AMD, built with ROCm 7.2 amdflang (validated)

Why: the cluster's ROCm-aware OpenMPI was built with gfortran, whose
`mpi_f08` module files cannot be used from amdflang. This build uses the
ROCm UCX 1.20.1 already installed. Validated with the mini-app on one and
two MI210s (2 ranks), including GPU-aware MPI with `rocm_ipc`.

```bash
R=/opt/ohpc/pub/utils/rocm/7.2/rocm-7.2.3
cd $HOME/sw/src
curl -sL -o openmpi-5.0.7.tar.bz2 \
  https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.7.tar.bz2
tar xjf openmpi-5.0.7.tar.bz2
cd openmpi-5.0.7
# clang 22 rejects an implicit declaration of memcpy in the ROCm component
sed -i '0,/#include "opal_config.h"/s//#include "opal_config.h"\n#include <string.h>/' \
  opal/mca/accelerator/rocm/accelerator_rocm_module.c
mkdir build-amdflang && cd build-amdflang
module purge          # no gnu15/openmpi5/libfabric in the environment
env -u LD_LIBRARY_PATH -u LIBRARY_PATH -u CPATH \
CC=$R/lib/llvm/bin/amdclang CXX=$R/lib/llvm/bin/amdclang++ FC=$R/lib/llvm/bin/amdflang \
CFLAGS='-Wno-error=int-conversion -Wno-error=implicit-function-declaration -Wno-error=incompatible-pointer-types' \
../configure --prefix=$HOME/sw/openmpi-5.0.7-rocm-amdflang \
  --with-ucx=/opt/ohpc/pub/mpi/rocm/ucx/1.20.1 --with-rocm=$R \
  --with-pmix --with-slurm --enable-mpi-fortran --disable-static \
  --disable-oshmem --with-hwloc=/opt/ohpc/pub/libs/hwloc \
  --without-ofi --without-verbs --without-psm2 --without-cuda
make -j 24 && make install
```

Notes for the module: it needs `hwloc/2.12.0` and `rocm/7.2` loaded (the
wrapper links `libhwloc.so.15`). The ROCm UCX 1.20.1 it uses has the
`rocm_copy`, `rocm_ipc`, shared memory and `tcp` transports but no
InfiniBand ones (`rc`, `ud`), so between nodes the GPU buffers go over
TCP: GHOST on 4 GPUs across a1 and a2 spends 95% of the step in the
exchange. For a system-wide install, build UCX with both `--with-rocm`
and `--with-verbs` (the gnu15 UCX 1.18 shows `rc_verbs`/`ud_verbs`, so
the nodes have the hardware), and use
`UCX_TLS=rocm_copy,rocm_ipc,sm,self,rc,tcp`. Recommended runtime environment for
device buffers: `OMPI_MCA_pml=ucx`, `UCX_TLS=rocm_copy,rocm_ipc,sm,self`.
Launch with `srun --mpi=pmix`.

## 2. LLVM 23.1.0 with flang and NVPTX OpenMP offload (validated)

Why: nvfortran 24.5 refuses OpenMP offload below compute capability 7.0
(Tesla P4 is 6.1), the ROCm LLVM has no NVPTX backend and gfortran 15 has
no nvptx offload. Validated with the mini-app on a P4 (1 and 2 ranks):
all checks pass. The build script is
`$HOME/sw/src/llvm-build/build_llvm.sbatch`; its content, minus the SLURM
header:

```bash
module purge
module load gnu15 cmake/4.1.2 cuda/12.4     # plus the default python module
SRC=$HOME/sw/src; VER=23.1.0; PREFIX=$HOME/sw/llvm-$VER
GCCDIR=/opt/ohpc/pub/compiler/gcc/15.2.0/lib/gcc/x86_64-pc-linux-gnu/15.2.0
cd $SRC
curl -sL -o llvm-project-$VER.src.tar.xz \
  https://github.com/llvm/llvm-project/releases/download/llvmorg-$VER/llvm-project-$VER.src.tar.xz
tar xJf llvm-project-$VER.src.tar.xz
mkdir -p llvm-build/build && cd llvm-build/build
# The freshly built clang must use GCC 15, not the devtoolset-7 it finds
# under /opt/rh, when it compiles the runtimes: config files next to the
# binaries, in the build tree and later in the install prefix.
mkdir -p bin
for c in clang clang++ flang clang-cpp; do
  echo "--gcc-install-dir=$GCCDIR" > bin/$c.cfg
done
cmake -G "Unix Makefiles" $SRC/llvm-project-$VER.src/llvm \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DLLVM_ENABLE_PROJECTS="clang;flang;lld;mlir" \
  -DLLVM_ENABLE_RUNTIMES="openmp;offload" \
  -DLLVM_TARGETS_TO_BUILD="X86;NVPTX" \
  -DLLVM_ENABLE_ASSERTIONS=OFF \
  -DLLVM_INCLUDE_TESTS=OFF -DLLVM_INCLUDE_EXAMPLES=OFF \
  -DLLVM_INCLUDE_BENCHMARKS=OFF -DLLVM_INCLUDE_DOCS=OFF \
  -DFLANG_INCLUDE_TESTS=OFF -DMLIR_INCLUDE_TESTS=OFF \
  -DLLVM_PARALLEL_LINK_JOBS=4 \
  -DLIBOMP_OMPD_SUPPORT=OFF -DLIBOMP_OMPD_GDB_SUPPORT=OFF \
  -DLIBOMPTARGET_DEVICE_ARCHITECTURES="sm_61;sm_70" \
  -DLLVM_RUNTIME_TARGETS="default;nvptx64-nvidia-cuda" \
  -DRUNTIMES_nvptx64-nvidia-cuda_LLVM_ENABLE_RUNTIMES="openmp" \
  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++
make -j 24 && make install
for c in clang clang++ flang clang-cpp; do
  echo "--gcc-install-dir=$GCCDIR" > $PREFIX/bin/$c.cfg
done
```

Notes: the GPU device runtime (`lib/nvptx64-nvidia-cuda/libompdevice.a`)
only exists when the NVPTX triple is in `LLVM_RUNTIME_TARGETS`. Adding
`libc` to the GPU runtimes needs PyYAML in the Python that CMake picks
(the default python module has it). Compiling for the P4 needs
`--offload-arch=sm_61` and, for `use_device_addr`,
`-fopenmp-version=52`; `ptxas` from `cuda/12.4` must be on PATH.

## 3. OpenMPI 5.0.7 with CUDA, built with the LLVM flang above (validated)

Built without UCX because the HPC-X UCX shipped with NVHPC has a
pkg-config file pointing at its build machine. Intra-node device buffers
then go through OpenMPI's own CUDA support.

```bash
L=$HOME/sw/llvm-23.1.0
C=/opt/ohpc/pub/compiler/nvhpc/24.5/Linux_x86_64/24.5/cuda/12.4   # or the cuda/12.4 module prefix
cd $HOME/sw/src/openmpi-5.0.7 && mkdir build-llvmflang-cuda && cd build-llvmflang-cuda
module purge; module load gnu15
env -u LIBRARY_PATH -u CPATH \
CC=$L/bin/clang CXX=$L/bin/clang++ FC=$L/bin/flang \
CFLAGS='-Wno-error=int-conversion -Wno-error=implicit-function-declaration -Wno-error=incompatible-pointer-types' \
../configure --prefix=$HOME/sw/openmpi-5.0.7-cuda-llvmflang \
  --with-cuda=$C --with-cuda-libdir=$C/lib64/stubs --without-ucx \
  --with-pmix --with-slurm --enable-mpi-fortran --disable-static \
  --disable-oshmem --with-hwloc=/opt/ohpc/pub/libs/hwloc \
  --without-ofi --without-verbs --without-psm2 --without-rocm
make -j 16 && make install
```

## 4. OSU micro-benchmarks 7.5 with ROCm (used for the bandwidth baseline)

```bash
cd $HOME/sw/src
curl -sL -o osu-micro-benchmarks-7.5.tar.gz \
  https://mvapich.cse.ohio-state.edu/download/mvapich/osu-micro-benchmarks-7.5.tar.gz
tar xzf osu-micro-benchmarks-7.5.tar.gz && cd osu-micro-benchmarks-7.5
module purge; module load rocm/7.2 hwloc/2.12.0
export PATH=$HOME/sw/openmpi-5.0.7-rocm-amdflang/bin:$PATH
./configure CC=mpicc CXX=mpicxx --prefix=$HOME/sw/osu-7.5-rocm --enable-rocm \
  --with-rocm=/opt/ohpc/pub/utils/rocm/7.2/rocm-7.2.3
make -j 8 -C c/mpi/pt2pt
```
