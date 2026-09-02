# -------------------------------------------------------------------
# GHOST CMakeLists
#
# If libraries are not in the path, use -DCMAKE_PREFIX_PATH
# and -DFFTW3_ROOT to set paths, e.g., from the build directory do
# "cmake .. -DCMAKE_PREFIX_PATH=/opt/openmpi-5.0.10 -DFFTW3_ROOT=/opt/fftw-3.3.10"
# Compilers can be declared in bash using "export CC=gcc;export FC=gfortran"
# You can define an install dir using -DCMAKE_INSTALL_PREFIX=/foo
# -DBIN_DIR=/foo/bin and -DLIB_DIR=/foo/lib place the executables and the
# Python library elsewhere than GHOST/bin and GHOST/lib (useful to keep
# several builds).
# Other Options:
# -DPRECISION=SINGLE(default)/DOUBLE
# -DP_HYBRID=ON (multithreaded, default: OFF)
# -DP_GPU=ON -DGPU_VENDOR=NVIDIA|AMD -DGPU_ARCH=<arch> (OpenMP offload,
#    default: OFF; needs nvfortran or flang, see src/README)
# -------------------------------------------------------------------

cmake_minimum_required(VERSION 3.12)
project(GHOST LANGUAGES Fortran C)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(UserConfig.cmake)

# -------------------------- Directories -----------------------------
# All Fortran sources go through the preprocessor: the offload build
# selects host or device loops with #if defined(GHOST_GPU) blocks.
set(CMAKE_Fortran_PREPROCESS ON)

set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(BIN_DIR "${CMAKE_SOURCE_DIR}/bin" CACHE PATH "Directory for the executables")
set(LIB_DIR "${CMAKE_SOURCE_DIR}/lib" CACHE PATH "Directory for the Python library")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}" ${CMAKE_MODULE_PATH})
file(MAKE_DIRECTORY "${BIN_DIR}")

# ----------------------- Compiler options ---------------------------
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O4 -funroll-loops")
    # GNU debug flags
    #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -gdwarf-2 -O0 \
    #      -fno-omit-frame-pointer -fbacktrace -fcheck=bounds \
    #      -finit-real=snan -finit-integer=-999999 -finit-logical=true \
    #      -finit-character=0 -Wall -Wextra")
    # Apply version-specific flags
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 \
           -fp-model fast=2 -align -heap-arrays")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
    # No LTO: it needs the gold plugin or lld and cannot build the shared
    # Python library with some flang versions, for little gain here
    # -fPIC explicitly: CMake does not apply POSITION_INDEPENDENT_CODE
    # to flang, and the Python library is shared
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -funroll-loops -fPIC \
           -Wno-incompatible-implicit-interfaces")
endif()

# ----------------------- Hybrid OpenMP ------------------------------
# P_HYBRID enables the hybrid MPI/OpenMP build. Besides the threaded
# loops it is what activates every '!$' conditional line, and so what
# makes main.fpp set nth from omp_get_max_threads() and fftp3D call
# fftp3d_init_threads and plan_with_nthreads: without it the FFTs run
# on a single thread. The flags are set globally, like the per-compiler
# options above, so that every library in the build tree gets them.
if(P_GPU AND NOT P_HYBRID)
    set(P_HYBRID ON CACHE BOOL "Enable hybrid OpenMP" FORCE)
    message(STATUS "P_GPU implies P_HYBRID: hybrid OpenMP enabled")
endif()
if(P_HYBRID)
    find_package(OpenMP COMPONENTS Fortran C)
    if(NOT OpenMP_Fortran_FOUND)
        # CMake's detection test fails with some flang versions although
        # the compiler supports OpenMP with the usual flag
        if(CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
            set(OpenMP_Fortran_FLAGS "-fopenmp")
            message(STATUS "OpenMP detection failed, assuming -fopenmp for flang")
        else()
            message(FATAL_ERROR "OpenMP is required for P_HYBRID and was not found")
        endif()
    endif()
    if(NOT OpenMP_C_FOUND AND CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C_FLAGS "-fopenmp")
    endif()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set(CMAKE_C_FLAGS       "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    message(STATUS "Hybrid OpenMP enabled (${OpenMP_Fortran_FLAGS})")
else()
    message(STATUS "Hybrid OpenMP disabled; use -DP_HYBRID=ON to enable")
endif()

# ----------------------- GPU offload --------------------------------
# P_GPU compiles the field kernels as OpenMP target regions and keeps
# the field arrays resident on the device (see src/memmgt/gdevice_mod
# and src/README). The macros GHOST_GPU and GPU_<VENDOR> select the
# device code paths and the vendor runtime entry points. Every library
# gets the flags because kernels are spread over several of them.
add_library(project_gpu INTERFACE)
if(P_GPU)
    if(NOT GPU_VENDOR MATCHES "^(NVIDIA|AMD)$")
        message(FATAL_ERROR "GPU_VENDOR must be NVIDIA or AMD (got '${GPU_VENDOR}')")
    endif()
    if(GPU_ARCH STREQUAL "")
        message(FATAL_ERROR "P_GPU needs GPU_ARCH (e.g. cc80 or sm_80 for NVIDIA, gfx90a for AMD)")
    endif()
    add_compile_definitions(GHOST_GPU GPU_${GPU_VENDOR})
    if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mp=gpu -gpu=${GPU_ARCH}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mp=gpu -gpu=${GPU_ARCH}")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -mp=gpu -gpu=${GPU_ARCH}")
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
        # flang checks directives against OpenMP 3.1 unless told otherwise
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} --offload-arch=${GPU_ARCH} -fopenmp-version=52")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp --offload-arch=${GPU_ARCH}")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fopenmp --offload-arch=${GPU_ARCH}")
    else()
        message(FATAL_ERROR "P_GPU needs nvfortran or flang (got ${CMAKE_Fortran_COMPILER_ID})")
    endif()
    # Vendor runtime, used by the allocator (hipMalloc/cudaMalloc) and
    # for device synchronization. FFT libraries are added by the FFTP
    # backend.
    if(GPU_VENDOR STREQUAL "AMD")
        set(ROCM_PATH "$ENV{ROCM_PATH}" CACHE PATH "ROCm install prefix")
        find_library(GPU_RUNTIME_LIB NAMES amdhip64 HINTS ${ROCM_PATH}/lib REQUIRED)
    else()
        set(CUDA_PATH "$ENV{CUDA_HOME}" CACHE PATH "CUDA install prefix")
        find_library(GPU_RUNTIME_LIB NAMES cudart HINTS ${CUDA_PATH}/lib64
                     ${CUDA_PATH}/lib64/stubs REQUIRED)
    endif()
    target_link_libraries(project_gpu INTERFACE ${GPU_RUNTIME_LIB})
    message(STATUS "GPU offload enabled: ${GPU_VENDOR} ${GPU_ARCH}, runtime ${GPU_RUNTIME_LIB}")
else()
    message(STATUS "GPU offload disabled; use -DP_GPU=ON to enable")
endif()

# ---------------------- Preprocessor definitions --------------------
# The .F90 sources are preprocessed by the Fortran compiler itself
# (CMAKE_Fortran_PREPROCESS is on for all sources); these definitions
# select the precision and the parameters of the parallel FFT.
add_compile_definitions(CSIZE_=${CSIZE} NSTRIP_=${NSTRIP} G${PRECISION}_PRECISION)

# ---------------------- MPI (Fortran) -------------------------------
# When the Fortran compiler is an MPI wrapper (mpifort, mpif90) MPI
# comes with it and CMake's detection is not needed (it fails with
# some compilers, e.g. flang with offload flags).
add_library(project_mpi INTERFACE)
if(CMAKE_Fortran_COMPILER MATCHES "mpif")
    message(STATUS "Fortran compiler is an MPI wrapper: ${CMAKE_Fortran_COMPILER}")
else()
    find_package(MPI REQUIRED Fortran)
    target_include_directories(project_mpi INTERFACE ${MPI_Fortran_INCLUDE_DIRS})
    target_link_libraries(project_mpi INTERFACE ${MPI_Fortran_LIBRARIES})
endif()

# ---------------------- FFT libraries -------------------------------
# fftp-3 uses FFTW on the host; fftp-gpu adds cuFFT or hipFFT on the
# device and keeps FFTW for the host path, so both need FFTW. With
# P_GPU the default backend is fftp-gpu.
add_library(project_fft INTERFACE)
if(P_GPU AND FFTP STREQUAL "fftp-3")
    set(FFTP "fftp-gpu" CACHE STRING "FFT library (fftp-3, fftp-gpu, fftp-mkl)" FORCE)
    message(STATUS "P_GPU: using the fftp-gpu backend")
endif()
if(FFTP STREQUAL "fftp-gpu" AND NOT P_GPU)
    message(FATAL_ERROR "FFTP=fftp-gpu requires -DP_GPU=ON")
endif()
if(FFTP STREQUAL "fftp-3" OR FFTP STREQUAL "fftp-gpu")
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")
    find_package(FFTW3 REQUIRED)
    target_include_directories(project_fft INTERFACE ${FFTW3_INCLUDE_DIRS})
    # With P_HYBRID, link FFTW's threading support. The _omp variant is
    # used rather than _threads so that FFTW shares the OpenMP thread
    # pool instead of starting pthreads of its own. It has to be listed
    # before the base library for static linking to resolve.
    if(P_HYBRID)
        if(PRECISION STREQUAL "SINGLE")
            set(_fftw_omp_name fftw3f_omp)
        else()
            set(_fftw_omp_name fftw3_omp)
        endif()
        find_library(FFTW3_OMP_LIBRARY NAMES ${_fftw_omp_name}
                     HINTS ${FFTW3_ROOT}/lib ${FFTW3_LIBRARY_DIRS}
                     PATH_SUFFIXES lib lib64 REQUIRED)
        target_link_libraries(project_fft INTERFACE ${FFTW3_OMP_LIBRARY})
        message(STATUS "FFTW threading: ${FFTW3_OMP_LIBRARY}")
    endif()
    if(PRECISION STREQUAL "SINGLE")
        target_link_libraries(project_fft INTERFACE ${FFTW3F_LIBRARIES} m)
    else()
        target_link_libraries(project_fft INTERFACE ${FFTW3_LIBRARIES} m)
    endif()
    if(FFTP STREQUAL "fftp-gpu")
        if(GPU_VENDOR STREQUAL "AMD")
            find_library(GPU_FFT_LIB NAMES hipfft HINTS ${ROCM_PATH}/lib REQUIRED)
        else()
            # cuFFT is in the CUDA toolkit or, with the NVIDIA HPC SDK,
            # in its math_libs tree next to the cuda directory
            file(GLOB _nv_mathlibs "${CUDA_PATH}/../../math_libs/*/lib64")
            find_library(GPU_FFT_LIB NAMES cufft HINTS ${CUDA_PATH}/lib64
                         ${CUDA_PATH}/../../math_libs/lib64 ${_nv_mathlibs} REQUIRED)
        endif()
        target_link_libraries(project_fft INTERFACE ${GPU_FFT_LIB})
        message(STATUS "Device FFT: ${GPU_FFT_LIB}")
    endif()
elseif(FFTP STREQUAL "fftp-mkl")
    set(MKLDIR "" CACHE PATH "MKL root")
    target_include_directories(project_fft INTERFACE "${MKLDIR}/include")
    target_link_directories(project_fft INTERFACE "${MKLDIR}/lib/intel64")
    target_link_libraries(project_fft INTERFACE
        mkl_intel_lp64 mkl_core mkl_sequential
    )
else()
    message(FATAL_ERROR "Unsupported FFT: ${FFTP}")
endif()

# ---------------------- IO library ---------------------------------
add_library(project_io INTERFACE)
if(NOT IOLIB STREQUAL "mpiio")
    message(FATAL_ERROR "GHOST supports only mpiio")
endif()

# ---------------------- Compile definitions -------------------------
set(SOLVER "GHOST" CACHE STRING "Builds GHOST with all solvers")
add_library(project_defs INTERFACE)
string(TOUPPER "${SOLVER}" SOLVER_UP)
target_compile_definitions(project_defs INTERFACE "${SOLVER_UP}_SOL")

# ---------------------- Include directories -------------------------
add_library(project_includes INTERFACE)
target_include_directories(project_includes INTERFACE
    "${SRC_DIR}"
    "${SRC_DIR}/utils"
    "${SRC_DIR}/pseudo"
    "${SRC_DIR}/memmgt"
    "${SRC_DIR}/pdes"
    "${SRC_DIR}/${FFTP}"
    "${SRC_DIR}/${IOLIB}"
    "${CMAKE_BINARY_DIR}"              # Generated .f90
    "${CMAKE_BINARY_DIR}/src/utils"    # Generated modules
    "${CMAKE_BINARY_DIR}/src/mpiio"
    "${CMAKE_BINARY_DIR}/src/pseudo"
    "${CMAKE_BINARY_DIR}/src/memmgt"
    "${CMAKE_BINARY_DIR}/src/pdes"
    "${CMAKE_BINARY_DIR}/src/particles"
    "${CMAKE_BINARY_DIR}/src/ic"
    "${CMAKE_BINARY_DIR}/src/force"
    "${CMAKE_BINARY_DIR}/src/icpart"
    "${CMAKE_BINARY_DIR}/src/steppers"
    "${CMAKE_BINARY_DIR}/src/${FFTP}"
)

# ---------------------- Subdirectories ------------------------------
add_subdirectory("${SRC_DIR}/${FFTP}")
add_subdirectory("${SRC_DIR}/pseudo")
add_subdirectory("${SRC_DIR}/utils")
add_subdirectory("${SRC_DIR}/mpiio")
add_subdirectory("${SRC_DIR}/memmgt")
add_subdirectory("${SRC_DIR}/pdes")
add_subdirectory("${SRC_DIR}/particles")
add_subdirectory("${SRC_DIR}/ic")
add_subdirectory("${SRC_DIR}/force")
add_subdirectory("${SRC_DIR}/icpart")
add_subdirectory("${SRC_DIR}/steppers")
add_subdirectory("${SRC_DIR}/tools")

# ---------------------- Main executable ------------------------------
add_executable(GHOST "${SRC_DIR}/main.F90")
target_link_libraries(GHOST
    PRIVATE
        project_includes
        project_gpu
        mpiio_lib
        particles_lib
        ic_lib
        force_lib
        icp_lib
        stepper_lib
)

set_target_properties(GHOST PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${BIN_DIR}"
)

# ---------------------- Installer ------------------------------------
include(GNUInstallDirs)
install(
    TARGETS GHOST BOOTS ghost_python
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# ---------------------------------------------------------------------
# End of CMakeLists.txt
# ---------------------------------------------------------------------
