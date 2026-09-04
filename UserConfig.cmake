# -------------------------------------------------------------------
# GHOST UserConfig CMake file
#
# If libraries are not in the path, use -DCMAKE_PREFIX_PATH
# and -DFFTW3_ROOT to set paths, e.g., from the build directory do
# "cmake .. -DCMAKE_PREFIX_PATH=/opt/openmpi-5.0.8 -DFFTW3_ROOT=/opt/fftw-3.3.10"
# -------------------------------------------------------------------

# ----------------------- Options / variants ------------------------
option(P_HYBRID "Enable hybrid OpenMP"                             OFF)
option(P_GPU    "Enable OpenMP offload to GPUs (implies P_HYBRID)" OFF)
set(GPU_VENDOR "NVIDIA" CACHE STRING "GPU vendor for P_GPU (NVIDIA or AMD)")
set(GPU_ARCH   ""       CACHE STRING "GPU architecture for P_GPU (e.g. cc80, sm_80, gfx90a)")
set(PRECISION  "SINGLE" CACHE STRING "Precision (SINGLE or DOUBLE)")
set(FFTP       "fftp-3" CACHE STRING "FFT library (fftp-3, fftp-gpu)")
set(IOLIB      "mpiio"  CACHE STRING "IO library (mpiio)")

# ----------------------- Advanced set up params --------------------
set(CSIZE  8 CACHE STRING "Cache size")
set(NSTRIP 1 CACHE STRING "Strip mining")
