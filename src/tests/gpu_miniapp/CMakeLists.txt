# -------------------------------------------------------------------
# GHOST GPU mini-app
#
# Standalone CMake project. Build with the MPI wrapper as the
# Fortran compiler, e.g.
#   FC=mpifort cmake -S . -B build -DGPU_VENDOR=NVIDIA -DGPU_ARCH=cc60 \
#        -DFFTW3_ROOT=/path/to/fftw
# Options:
#   GPU_VENDOR        NONE (host OpenMP), NVIDIA (nvfortran, cuFFT),
#                     AMD (amdflang, hipFFT)
#   GPU_ARCH          cc60, cc80, gfx90a, ...
#   PRECISION         SINGLE or DOUBLE
#   GPU_LOOP_SPELLING 1: target teams distribute parallel do
#                     2: target teams loop
#   TEST_DIRECT_DT    compile the kernel that indexes derived type
#                     components directly (ON)
#   USE_VENDOR_ALLOC  device copies from hipMalloc/cudaMalloc instead of
#                     omp_target_alloc (ON: memory from omp_target_alloc
#                     takes a path 30x slower through GPU-aware MPI)
#   FFTW3_ROOT        FFTW install prefix (host reference transform)
#   ROCM_PATH         ROCm prefix (AMD), default from the environment
#   CUDA_PATH         CUDA prefix (NVIDIA with LLVM flang)
# -------------------------------------------------------------------
cmake_minimum_required(VERSION 3.20)
project(gpu_miniapp LANGUAGES Fortran)

set(GPU_VENDOR "NONE" CACHE STRING "NONE, NVIDIA or AMD")
set(GPU_ARCH "" CACHE STRING "GPU architecture (cc60, gfx90a, ...)")
set(PRECISION "SINGLE" CACHE STRING "SINGLE or DOUBLE")
set(GPU_LOOP_SPELLING 1 CACHE STRING "1 distribute parallel do, 2 teams loop")
option(TEST_DIRECT_DT "Compile the direct derived-type kernel" ON)
option(USE_VENDOR_ALLOC "Allocate device copies with hipMalloc/cudaMalloc instead of omp_target_alloc" ON)
set(FFTW3_ROOT "" CACHE PATH "FFTW3 prefix")
set(ROCM_PATH "$ENV{ROCM_PATH}" CACHE PATH "ROCm prefix")
set(CUDA_PATH "$ENV{CUDA_HOME}" CACHE PATH "CUDA prefix (LLVM flang on NVIDIA)")

add_executable(gpu_miniapp mini_mod.F90 mini_kernels.F90 mini_fft.F90 mini_main.F90)
target_include_directories(gpu_miniapp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set_source_files_properties(mini_mod.F90 mini_kernels.F90 mini_fft.F90 mini_main.F90
    PROPERTIES Fortran_PREPROCESS ON)

# ---- definitions
target_compile_definitions(gpu_miniapp PRIVATE
    G${PRECISION}_PRECISION GPU_LOOP_SPELLING=${GPU_LOOP_SPELLING})
if(TEST_DIRECT_DT)
    target_compile_definitions(gpu_miniapp PRIVATE TEST_DIRECT_DT)
endif()
if(USE_VENDOR_ALLOC)
    target_compile_definitions(gpu_miniapp PRIVATE USE_VENDOR_ALLOC)
endif()
if(NOT GPU_VENDOR STREQUAL "NONE")
    target_compile_definitions(gpu_miniapp PRIVATE GHOST_GPU GPU_${GPU_VENDOR})
endif()

# ---- compiler flags
message(STATUS "Fortran compiler id: ${CMAKE_Fortran_COMPILER_ID}")
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(gpu_miniapp PRIVATE -O2 -fopenmp -ffree-line-length-none)
    target_link_options(gpu_miniapp PRIVATE -fopenmp)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    target_compile_options(gpu_miniapp PRIVATE -O2 -Minfo=mp)
    if(GPU_VENDOR STREQUAL "NVIDIA")
        target_compile_options(gpu_miniapp PRIVATE -mp=gpu -gpu=${GPU_ARCH})
        target_link_options(gpu_miniapp PRIVATE -mp=gpu -gpu=${GPU_ARCH} -cudalib=cufft)
    else()
        target_compile_options(gpu_miniapp PRIVATE -mp)
        target_link_options(gpu_miniapp PRIVATE -mp)
    endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
    # upstream flang checks directives against OpenMP 3.1 unless told
    # otherwise; use_device_addr needs 5.0 and later
    target_compile_options(gpu_miniapp PRIVATE -O2 -fopenmp -fopenmp-version=52)
    target_link_options(gpu_miniapp PRIVATE -fopenmp)
    if(GPU_VENDOR STREQUAL "AMD")
        target_compile_options(gpu_miniapp PRIVATE --offload-arch=${GPU_ARCH})
        target_link_options(gpu_miniapp PRIVATE --offload-arch=${GPU_ARCH})
        target_link_libraries(gpu_miniapp PRIVATE
            ${ROCM_PATH}/lib/libhipfft.so ${ROCM_PATH}/lib/libamdhip64.so)
    elseif(GPU_VENDOR STREQUAL "NVIDIA")
        # LLVM flang with the NVPTX backend (e.g. GPU_ARCH=sm_61)
        target_compile_options(gpu_miniapp PRIVATE --offload-arch=${GPU_ARCH})
        target_link_options(gpu_miniapp PRIVATE --offload-arch=${GPU_ARCH})
        find_library(CUFFT_LIB NAMES cufft HINTS ${CUDA_PATH}/lib64
                     ${CUDA_PATH}/../math_libs/lib64 REQUIRED)
        find_library(CUDART_LIB NAMES cudart HINTS ${CUDA_PATH}/lib64 REQUIRED)
        target_link_libraries(gpu_miniapp PRIVATE ${CUFFT_LIB} ${CUDART_LIB})
        message(STATUS "cuFFT: ${CUFFT_LIB}")
    endif()
else()
    message(FATAL_ERROR "Unsupported Fortran compiler ${CMAKE_Fortran_COMPILER_ID}")
endif()

# ---- FFTW (host reference; C library, usable from any compiler)
if(PRECISION STREQUAL "DOUBLE")
    set(_fftw_name fftw3)
else()
    set(_fftw_name fftw3f)
endif()
find_library(FFTW3_LIB NAMES ${_fftw_name} HINTS ${FFTW3_ROOT}/lib ${FFTW3_ROOT}/lib64 REQUIRED)
target_link_libraries(gpu_miniapp PRIVATE ${FFTW3_LIB} m)
message(STATUS "FFTW: ${FFTW3_LIB}")
