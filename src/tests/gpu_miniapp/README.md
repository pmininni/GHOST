# GHOST GPU mini-app

Standalone probe for the OpenMP offload design of GHOST (phase 1 of
the GPU plan). It reproduces, in about 1500 lines, the structures of
the code that matter for device residency and checks, per compiler
and vendor, the assumptions the design rests on:

1. **Dual-pragma kernels** on explicit-shape dummy arrays that are
   already resident on the device: `!$omp target teams distribute
   parallel do collapse(3)` (or `target teams loop`) in an offload
   build, `!$omp parallel do collapse(2)` over `DO CONCURRENT` in a
   host build. The loop headers live in `kloop3_begin.h` and
   `rloop3_begin.h`.
2. **Residency at the allocation choke points**: `GState_alloc`
   (states), `initialize_pool` (workspace) and `fftp3d_create_plan`
   (FFT buffers) do `target enter data map(alloc:)`; nothing else maps.
3. **Arrays of derived types with allocatable components** updated
   three ways: passing the component as an argument, through a
   pointer alias, and indexing `u(ic)%ccomp(k,j,i)` directly inside
   the kernel (the last one can be excluded with `-DTEST_DIRECT_DT=OFF`).
4. **Workspace pool pointers** resolving in the device table
   (`omp_target_is_present`).
5. **Reductions** on the device, compared with the host and checked
   for run-to-run bitwise equality.
6. **cuFFT / hipFFT** on resident arrays through ISO_C_BINDING and
   `use_device_addr`, compared with FFTW on the host; the full
   parallel 3D transform (2D FFT, slab exchange, transpose, 1D FFT).
7. **The slab exchange** with GPU-aware MPI, in three variants
   selectable at run time: MPI derived datatypes on device buffers,
   device-side pack into contiguous buffers, and host staging. All are
   checked and timed.

## Layout

| file | content |
|---|---|
| `gpu_defs.h` | macro header: FFT library names and handle type per vendor and precision |
| `kloop3_begin.h`, `rloop3_begin.h` | dual-pragma loop headers (spectral and real blocks) |
| `mini_mod.F90` | precision, MPI ranges, grid, wavenumbers, device helpers, state type, workspace pool |
| `mini_kernels.F90` | kernels and their serial host references |
| `mini_fft.F90` | plans, exchanges, transposes, parallel transforms |
| `mini_main.F90` | driver: checks and timings |
| `build_*.sh`, `run_*.sh` | per-vendor build and SLURM run scripts |

## Building and running

```
./build_host.sh                 # gfortran + OpenMPI, host OpenMP (harness check)
./build_nvidia.sh               # nvfortran + HPC-X (needs a >= cc70 GPU, see below)
./build_amd.sh                  # amdflang + ROCm-aware OpenMPI built with amdflang
NP=2 NG=1 ./run_amd.sh 128 10   # N=128, 10 timing iterations, 2 ranks on 1 GPU
```

Extra CMake options go to the build scripts, e.g.
`./build_amd.sh -DGPU_LOOP_SPELLING=2 -DPRECISION=DOUBLE`.

To see whether any data moves during the kernels, run with
`LIBOMPTARGET_INFO=-1` (LLVM/AMD) or `NVCOMPILER_ACC_NOTIFY=3`
(NVHPC): only the explicit `target update` calls of the driver should
appear.

## Findings on the sakura cluster

### Toolchains

- **nvfortran 24.5 does not offload with OpenMP below compute
  capability 7.0** (`NVFORTRAN-W-0155`), with either the CUDA 12.4 or
  the CUDA 11.8 toolchain. The Tesla P4 (cc 6.1) in g1 and g2 therefore
  cannot run the OpenMP path with the NVIDIA compiler; OpenACC does
  compile for cc60. nvfortran compiles the whole mini-app for cc70
  without a single diagnostic. The ROCm LLVM has only the AMD backend
  and gfortran 15 has no nvptx offload, so an upstream LLVM flang with
  the NVPTX backend (`$HOME/sw/llvm-23.1.0`) is the only OpenMP route to
  those cards.
- The ROCm-aware OpenMPI installed on the cluster was built with
  gfortran, so its `mpi_f08` module cannot be used from amdflang. A
  second OpenMPI 5.0.7, built with amdflang against the ROCm UCX 1.20.1,
  is used instead (`$HOME/sw/openmpi-5.0.7-rocm-amdflang`; needs
  `--disable-oshmem` and a `<string.h>` include in the ROCm accelerator
  component with clang 22). Jobs must be launched with `srun
  --mpi=pmix`, otherwise every task becomes its own MPI world.
- Host build (gfortran, 2 ranks): all checks pass, both exchange
  variants give bitwise identical transforms.

### AMD MI210, ROCm 7.2 amdflang (flang 22)

All kernels, the argument and pointer forms of the state update, the
device reduction (bitwise reproducible run to run), and the hipFFT
transforms pass, in single and double precision and with both loop
spellings. Three things had to be learned the hard way:

1. **`target enter data` on a derived-type component maps the chain of
   enclosing descriptors.** For `this%real_entries_(i)%array` inside a
   type-bound procedure, flang also mapped the 40-byte CLASS descriptor
   of `this`, which lives on the caller's stack. That entry stayed in
   the device table after the return and later collided with the
   transient descriptors kernels map ("explicit extension not
   allowed"). The mini-app therefore creates device copies at the three
   allocation choke points with `omp_target_alloc` and
   `omp_target_associate_ptr` on the raw data address
   (`dev_alloc_c/r` in `mini_mod.F90`), so only data ever sits in the
   table. Kernels then find their explicit-shape dummies present and
   transfer nothing.
2. **Several target regions expanded from the same `#include` line in
   one routine get their device kernels swapped.** flang names offload
   entries by source location (`derivk3_l8`, `derivk3_l8_1`,
   `derivk3_l8_2`); with the loop header in an include file the three
   `derivk3` kernels were launched in the wrong order and the x and z
   derivatives came out as z and x. The directives are now written out
   in every kernel.
3. **Indexing `u(ic)%ccomp(k,j,i)` directly inside a kernel silently
   does nothing to the resident data.** flang maps the array of
   `GStateComp` by copying its three 96-byte descriptors to the device
   verbatim, with base addresses still pointing at host memory. The
   result on the device is unchanged and no error is raised. Passing the
   component as an argument or through a pointer alias both work, so the
   loops in `canuto_stepper.f90` and at the end of `dudt` that index
   components directly must be rewritten one of those two ways.

4. **Memory from `omp_target_alloc` travels slowly through GPU-aware
   MPI**, and a rank sending a device buffer to itself through MPI is
   slower still. Device copies are therefore allocated with `hipMalloc`
   (or `cudaMalloc`) and associated with the host arrays
   (`USE_VENDOR_ALLOC`, on by default), and the diagonal block of the
   transpose is copied with a device kernel instead of MPI. The vendor
   runtime's current device must also be set to the OpenMP device
   (`hipSetDevice`), or hipFFT plans and allocations land on device 0
   for every rank and the two-GPU run faults.
5. **MPI derived datatypes on device buffers are unusable** (seconds
   per exchange at N=256). The exchange has to pack on the device into
   contiguous messages. The right transport is UCX with the ROCm
   transports (`OMPI_MCA_pml=ucx`, `UCX_TLS=rocm_copy,rocm_ipc,sm,self`);
   the default ob1 path is about 10 times slower for device buffers. The
   OSU benchmark with this MPI reaches 75 GB/s device to device between
   two ranks with UCX, 1.2 GB/s with ob1.

### Performance, AMD MI210

Milliseconds per forward+inverse transform pair, 2 ranks, 2 GPUs, one
host core per rank, single precision, `hipMalloc` copies, UCX ROCm
transports. "host" is the FFTW path on one core per rank.

| N   | path                         | total  | FFTs  | exchange | transposes |
|-----|------------------------------|-------:|------:|---------:|-----------:|
| 128 | host (FFTW, datatypes)       |  47.5  | 30.4  |   8.4    |   9.8      |
| 128 | device, datatypes            | 783    |  0.39 | 782      |   0.25     |
| 128 | device, packed               |   0.88 |  0.19 |   0.54   |   0.16     |
| 128 | device, host staged          |   9.3  |  0.30 |   8.8    |   0.22     |
| 256 | host (FFTW, datatypes)       | 415    | 271   |  56      |  93        |
| 256 | device, datatypes            | 3737   |  0.81 | 3735     |   1.4      |
| 256 | device, packed               |   3.83 |  0.63 |   1.87   |   1.35     |
| 256 | device, host staged          |  72    |  0.82 |  70      |   1.4      |

The three spectral kernels (derivative, Laplacian, dealiased update)
take 0.6 ms together at N=256. Device reductions are bitwise
reproducible from run to run. With `omp_target_alloc` memory the packed
exchange at N=128 takes 28 ms instead of 0.5 ms.

### NVIDIA Tesla P4 with upstream LLVM flang 23.1

Toolchain: `$HOME/sw/llvm-23.1.0` (clang, flang, lld, mlir, the OpenMP
offload runtime, and the NVPTX device runtime obtained by adding
`nvptx64-nvidia-cuda` to `LLVM_RUNTIME_TARGETS`), with a CUDA-aware
OpenMPI 5.0.7 built with that flang (`$HOME/sw/openmpi-5.0.7-cuda-llvmflang`,
without UCX: the HPC-X UCX pkg-config file points at build-machine paths).
Two source-level accommodations were needed to compile:

- flang checks directives against OpenMP 3.1 by default, so
  `use_device_addr` needs `-fopenmp-version=52` (or 50 and later).
- `ABS` of a complex inside a target region is lowered to a call with a
  complex argument, which flang 23.1 cannot yet do for NVPTX ("not yet
  implemented: handle complex argument types"). The squared modulus is
  written as `REAL(z)**2+AIMAG(z)**2`. GHOST's `energy` and spectra
  routines use `abs(z)**2` and will need the same rewrite.

**Results (job 5777, g2):** every check passes with 1 and 2 ranks at
N=64 and N=128: kernels, pool pointers, argument and pointer updates,
bitwise-reproducible reductions, cuFFT through `use_device_addr`, and all
three exchange variants (the direct derived-type kernel is excluded with
`-DTEST_DIRECT_DT=OFF`; with it compiled in, the copied descriptors make
the kernel fault with an illegal memory access instead of silently doing
nothing as on AMD). Milliseconds per transform pair, 2 ranks sharing the
one P4 allocated to the job, N=128:

| path                    | total | FFTs | exchange | transposes |
|-------------------------|------:|-----:|---------:|-----------:|
| host (FFTW, datatypes)  | 35.8  | 27.9 |   4.4    |   4.9      |
| device, datatypes       | 2228  |  2.2 | 2216     |  10.9      |
| device, packed          | 27.9  |  2.4 |  14.7    |  11.5      |
| device, host staged     | 20.6  |  1.6 |   8.2    |  11.8      |

cuFFT is 12 times faster than FFTW on one core, as expected, but the
OpenMP kernels generated by flang 23.1 for sm_61 are slow: the three
spectral kernels take 25 ms (0.6 ms on the MI210) and the transposes 11
ms. The launch geometry is not the reason: every kernel runs as 3200
blocks of 128 threads in SPMD mode, and the `teams loop` spelling gives
the same times. `OMP_NUM_TEAMS=80 OMP_TEAMS_THREAD_LIMIT=256` halves the
kernel time (12.4 ms), which still leaves a factor of 20 to the MI210
against a hardware ratio of about 8 in memory bandwidth. What remains is
the quality of flang's NVPTX code for descriptor-indexed collapsed loops,
a matter for the production compiler on production hardware, not for the
design. The P4 path is kept as a correctness check of the NVIDIA build.

### Not yet available

- Node a2 lacked the system `libdrm` libraries until 2026-09-02; it is
  now identical to a1.
