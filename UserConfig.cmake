# -------------------------------------------------------------------
# GHOST UserConfig CMake file
#
# If libraries are not in the path, use -DCMAKE_PREFIX_PATH
# and -DFFTW3_ROOT to set paths, e.g., from the build directory do
# "cmake .. -DCMAKE_PREFIX_PATH=/opt/openmpi-5.0.8 -DFFTW3_ROOT=/opt/fftw-3.3.10"
# -------------------------------------------------------------------

# ----------------------- Advanced set up params --------------------
set(IKIND 8  CACHE STRING "Integer kind")
set(CSIZE 8  CACHE STRING "Cache size")
set(NSTRIP 1 CACHE STRING "Strip mining")

# ----------------------- Options / variants ------------------------
option(ARBSIZE  "Enable arbitrary box sizes" OFF)
option(P_HYBRID "Enable hybrid OpenMP"       OFF)
set(SOLVER "HD"        CACHE STRING "Solver type (HD, MHD)")
set(PRECISION "SINGLE" CACHE STRING "Precision (SINGLE or DOUBLE)")
set(FFTP "fftp-3"      CACHE STRING "FFT library (fftp-3, fftp-mkl)")
set(IOLIB "mpiio"      CACHE STRING "IO library (mpiio)")
