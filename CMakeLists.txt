# -------------------------------------------------------------------
# GHOST CMakeLists
#
# If libraries are not in the path, use -DCMAKE_PREFIX_PATH
# and -DFFTW3_ROOT to set paths, e.g., from the build directory do
# "cmake .. -DCMAKE_PREFIX_PATH=/opt/openmpi-5.0.10 -DFFTW3_ROOT=/opt/fftw-3.3.10"
# Compilers can be declared in bash using "export CC=gcc;export FC=gfortran"
# You can define an install dir using -DCMAKE_INSTALL_PREFIX=/foo
# Other Options:
# -DPRECISION=SINGLE(default)/DOUBLE
# -DP_HYBRID=ON (multithreaded, default: OFF)
# -------------------------------------------------------------------

cmake_minimum_required(VERSION 3.12)
project(GHOST LANGUAGES Fortran C)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(UserConfig.cmake)

# -------------------------- Directories -----------------------------
set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(BIN_DIR "${CMAKE_SOURCE_DIR}/bin")
set(LIB_DIR "${CMAKE_SOURCE_DIR}/lib")
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
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -flto -funroll-loops \
           -Wno-incompatible-implicit-interfaces")
    endif()

# ----------------------- Hybrid OpenMP ------------------------------
# P_HYBRID enables the hybrid MPI/OpenMP build. Besides the threaded
# loops it is what activates every '!$' conditional line, and so what
# makes main.fpp set nth from omp_get_max_threads() and fftp3D call
# fftp3d_init_threads and plan_with_nthreads: without it the FFTs run
# on a single thread. The flags are set globally, like the per-compiler
# options above, so that every library in the build tree gets them.
if(P_HYBRID)
    find_package(OpenMP REQUIRED COMPONENTS Fortran C)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set(CMAKE_C_FLAGS       "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    message(STATUS "Hybrid OpenMP enabled (${OpenMP_Fortran_FLAGS})")
else()
    message(STATUS "Hybrid OpenMP disabled; use -DP_HYBRID=ON to enable")
endif()

# ---------------------- Preprocessing function ----------------------
# Preprocess .fpp → .f90 using the C preprocessor (in the build tree)
function(preprocess_fpp OUTVAR FPP_PATH)
    get_filename_component(fname "${FPP_PATH}" NAME_WE)
    set(out "${CMAKE_BINARY_DIR}/${fname}.f90")
    add_custom_command(
        OUTPUT "${out}"
        COMMAND ${CMAKE_C_COMPILER} -E
            -DCSIZE_=${CSIZE} -DNSTRIP_=${NSTRIP}
            -DG${PRECISION}_PRECISION "${FPP_PATH}" -o "${out}"
        DEPENDS "${FPP_PATH}"
        COMMENT "Preprocessing ${fname}.fpp → build/${fname}.f90"
        VERBATIM
    )
    set(${OUTVAR} "${out}" PARENT_SCOPE)
endfunction()

# ---------------------- MPI (Fortran) -------------------------------
find_package(MPI REQUIRED Fortran)
add_library(project_mpi INTERFACE)
target_include_directories(project_mpi INTERFACE ${MPI_Fortran_INCLUDE_DIRS})
target_link_libraries(project_mpi INTERFACE ${MPI_Fortran_LIBRARIES})

# ---------------------- FFT libraries -------------------------------
add_library(project_fft INTERFACE)
if(FFTP STREQUAL "fftp-3")
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

# ---------------------- Preprocess FPP sources ----------------------
preprocess_fpp(GUTILS_SRC   "${SRC_DIR}/utils/gutils_mod.fpp")
preprocess_fpp(SPLINE_SRC   "${SRC_DIR}/particles/gpsplineint_mod.fpp")
preprocess_fpp(PRECIS_SRC   "${SRC_DIR}/${FFTP}/fprecision_mod.fpp")
preprocess_fpp(FFTMOD_SRC   "${SRC_DIR}/${FFTP}/fftp_mod.fpp")
preprocess_fpp(FFTP3D_SRC   "${SRC_DIR}/${FFTP}/fftp3D.fpp")
preprocess_fpp(MAIN_SRC     "${SRC_DIR}/main.fpp")

add_custom_target(project_preprocess
    DEPENDS "${GUTILS_SRC}" "${SPLINE_SRC}" "${PRECIS_SRC}"
    "${FFTMOD_SRC}" "${FFTP3D_SRC}"
)   

# ---------------------- Main executable ------------------------------
add_executable(GHOST ${MAIN_SRC})
add_dependencies(GHOST project_preprocess)
target_link_libraries(GHOST
    PRIVATE
        project_includes
        mpiio_lib
        memmgt_lib
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
