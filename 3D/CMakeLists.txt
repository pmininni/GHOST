# -------------------------------------------------------------------
# GHOST CMakeLists
#
# If libraries are not in the path, use -DCMAKE_PREFIX_PATH
# and -DFFTW3_ROOT to set paths, e.g., from the build directory do
# "cmake .. -DCMAKE_PREFIX_PATH=/opt/openmpi-5.0.8 -DFFTW3_ROOT=/opt/fftw-3.3.10"
# -------------------------------------------------------------------

cmake_minimum_required(VERSION 3.12)
project(GHOST LANGUAGES Fortran C)

include(UserConfig.cmake)

# -------------------------- Directories -----------------------------
set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(BIN_DIR "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}" ${CMAKE_MODULE_PATH})
file(MAKE_DIRECTORY "${BIN_DIR}")

# ----------------------- Compiler options ---------------------------
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # Apply flags for all Fortran files
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -funroll-loops")
    # Apply version-specific flags
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "9.0")
    message(STATUS "GNU Fortran version is 9.0 or newer, enabling -fallow-argument-mismatch.")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w -fallow-argument-mismatch")
    endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -no-prec-div
           -fp-model fast=2 -align -heap-arrays"
    )
    endif()

# ---------------------- Preprocessing function ----------------------
# Preprocess .fpp → .f90 using the C preprocessor (in the build tree)
function(preprocess_fpp OUTVAR FPP_PATH)
    get_filename_component(fname "${FPP_PATH}" NAME_WE)
    set(out "${CMAKE_BINARY_DIR}/${fname}.f90")
    add_custom_command(
        OUTPUT "${out}"
        COMMAND ${CMAKE_C_COMPILER} -E
            -DNX_=${NX} -DNY_=${NY} -DNZ_=${NZ} -DORD_=${ORD}
            -DIKIND_=${IKIND} -DCSIZE_=${CSIZE} -DNSTRIP_=${NSTRIP}
            -DG${PRECISION}_PRECISION -D${SOLVER}_SOL
            "${FPP_PATH}" -o "${out}"
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
add_library(project_defs INTERFACE)

string(TOUPPER "${SOLVER}" SOLVER_UP)
target_compile_definitions(project_defs INTERFACE "${SOLVER_UP}_SOL")

if(ARBSIZE)
    target_compile_definitions(project_defs INTERFACE DEF_ARBSIZE_)
endif()

# ---------------------- Include directories -------------------------
add_library(project_includes INTERFACE)
target_include_directories(project_includes INTERFACE
    "${SRC_DIR}"
    "${SRC_DIR}/include"
    "${SRC_DIR}/utils"
    "${SRC_DIR}/pseudo"
    "${SRC_DIR}/${FFTP}"
    "${SRC_DIR}/${IOLIB}"
    "${CMAKE_BINARY_DIR}"              # Generated .f90
    "${CMAKE_BINARY_DIR}/src/utils"    # Generated modules
    "${CMAKE_BINARY_DIR}/src/mpiio"
    "${CMAKE_BINARY_DIR}/src/${FFTP}"
)

# ---------------------- Subdirectories ------------------------------
add_subdirectory("${SRC_DIR}/utils")
add_subdirectory("${SRC_DIR}/mpiio")
add_subdirectory("${SRC_DIR}/${FFTP}")

# ---------------------- Preprocess FPP sources ----------------------
preprocess_fpp(PSEUDO3D_SRC "${SRC_DIR}/pseudo/pseudospec3D_mod.fpp")
preprocess_fpp(GUTILS_SRC   "${SRC_DIR}/utils/gutils_mod.fpp")
preprocess_fpp(PRECIS_SRC   "${SRC_DIR}/${FFTP}/fprecision_mod.fpp")
preprocess_fpp(FFTMOD_SRC   "${SRC_DIR}/${FFTP}/fftp_mod.fpp")
preprocess_fpp(FFTP3D_SRC   "${SRC_DIR}/${FFTP}/fftp3D.fpp")
preprocess_fpp(MAIN_SRC     "${SRC_DIR}/main3D.fpp")

add_custom_target(project_preprocess
    DEPENDS "${PSEUDO3D_SRC}" "${GUTILS_SRC}" "${PRECIS_SRC}"
    "${FFTMOD_SRC}" "${FFTP3D_SRC}"
)   

# ---------------------- Source selection -----------------------------
set(HD_SRC
    "${PSEUDO3D_SRC}"
    "${SRC_DIR}/pseudo/module_dns.f90"
    "${SRC_DIR}/pseudo/pseudospec3D_hd.f90"
    "${MAIN_SRC}"
)

set(MHD_SRC
    ${HD_SRC}
    "${SRC_DIR}/pseudo/pseudospec3D_mhd.f90"
)

if(SOLVER STREQUAL "HD")
    set(SOLVER_SRC ${HD_SRC})
elseif(SOLVER STREQUAL "MHD")
    set(SOLVER_SRC ${MHD_SRC})
else()
    message(FATAL_ERROR "SOLVER must be HD or MHD")
endif()

# ---------------------- Main executable ------------------------------
add_library(pseudo_lib STATIC "${PSEUDO3D_SRC}")
add_dependencies(pseudo_lib fftp_lib)
target_include_directories(pseudo_lib PUBLIC
    ${CMAKE_BINARY_DIR}/src/${FFTP}
)
add_executable(${SOLVER} ${SOLVER_SRC})
add_dependencies(${SOLVER} project_preprocess timer_lib fftp_lib
    utils_lib pseudo_lib mpiio_lib
)
target_link_libraries(${SOLVER}
    PRIVATE
        project_mpi
        project_fft
        project_includes
	timer_lib
	fftp_lib
	pseudo_lib
	utils_lib
	mpiio_lib
)

set_target_properties(${SOLVER} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${BIN_DIR}"
)

# ---------------------------------------------------------------------
# End of CMakeLists.txt
# ---------------------------------------------------------------------
