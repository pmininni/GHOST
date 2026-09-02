!=================================================================
! gpu_defs.h
! Preprocessor definitions for the GHOST GPU mini-app.
!
! Build-time macros (set by CMake):
!   GHOST_GPU           : offload build (else host OpenMP build)
!   GPU_NVIDIA/GPU_AMD  : vendor (selects cuFFT/hipFFT and HIP/CUDA
!                         runtime entry points)
!   GSINGLE_PRECISION/GDOUBLE_PRECISION
!   GPU_LOOP_SPELLING   : 1 = target teams distribute parallel do
!                         2 = target teams loop
!
! Only string literals and type keywords are defined here, no token
! pasting, so that the gfortran, nvfortran and flang preprocessors
! all agree on the result.
!=================================================================

! ------------------------- FFTW (host reference) -----------------
#if defined(GDOUBLE_PRECISION)
#  define C_FFTW_PLAN_MANY_R2C "fftw_plan_many_dft_r2c"
#  define C_FFTW_PLAN_MANY_C2R "fftw_plan_many_dft_c2r"
#  define C_FFTW_PLAN_MANY_C2C "fftw_plan_many_dft"
#  define C_FFTW_EXEC_R2C      "fftw_execute_dft_r2c"
#  define C_FFTW_EXEC_C2R      "fftw_execute_dft_c2r"
#  define C_FFTW_EXEC_C2C      "fftw_execute_dft"
#  define C_FFTW_DESTROY       "fftw_destroy_plan"
#else
#  define C_FFTW_PLAN_MANY_R2C "fftwf_plan_many_dft_r2c"
#  define C_FFTW_PLAN_MANY_C2R "fftwf_plan_many_dft_c2r"
#  define C_FFTW_PLAN_MANY_C2C "fftwf_plan_many_dft"
#  define C_FFTW_EXEC_R2C      "fftwf_execute_dft_r2c"
#  define C_FFTW_EXEC_C2R      "fftwf_execute_dft_c2r"
#  define C_FFTW_EXEC_C2C      "fftwf_execute_dft"
#  define C_FFTW_DESTROY       "fftwf_destroy_plan"
#endif

! ------------------------- Device FFT ----------------------------
#if defined(GPU_NVIDIA)
#  define GFFT_HANDLE_T        INTEGER(C_INT)
#  define C_GFFT_PLAN_MANY     "cufftPlanMany"
#  define C_GFFT_DESTROY       "cufftDestroy"
#  define C_DEV_SYNC           "cudaDeviceSynchronize"
#  define C_DEV_MALLOC         "cudaMalloc"
#  define C_DEV_SETDEVICE      "cudaSetDevice"
#  define C_DEV_FREE           "cudaFree"
#  if defined(GDOUBLE_PRECISION)
#    define C_GFFT_EXEC_R2C    "cufftExecD2Z"
#    define C_GFFT_EXEC_C2R    "cufftExecZ2D"
#    define C_GFFT_EXEC_C2C    "cufftExecZ2Z"
#  else
#    define C_GFFT_EXEC_R2C    "cufftExecR2C"
#    define C_GFFT_EXEC_C2R    "cufftExecC2R"
#    define C_GFFT_EXEC_C2C    "cufftExecC2C"
#  endif
#elif defined(GPU_AMD)
#  define GFFT_HANDLE_T        TYPE(C_PTR)
#  define C_GFFT_PLAN_MANY     "hipfftPlanMany"
#  define C_GFFT_DESTROY       "hipfftDestroy"
#  define C_DEV_SYNC           "hipDeviceSynchronize"
#  define C_DEV_MALLOC         "hipMalloc"
#  define C_DEV_SETDEVICE      "hipSetDevice"
#  define C_DEV_FREE           "hipFree"
#  if defined(GDOUBLE_PRECISION)
#    define C_GFFT_EXEC_R2C    "hipfftExecD2Z"
#    define C_GFFT_EXEC_C2R    "hipfftExecZ2D"
#    define C_GFFT_EXEC_C2C    "hipfftExecZ2Z"
#  else
#    define C_GFFT_EXEC_R2C    "hipfftExecR2C"
#    define C_GFFT_EXEC_C2R    "hipfftExecC2R"
#    define C_GFFT_EXEC_C2C    "hipfftExecC2C"
#  endif
#endif

! cuFFT and hipFFT share the transform type codes
#if defined(GDOUBLE_PRECISION)
#  define GFFT_TYPE_R2C  106
#  define GFFT_TYPE_C2R  108
#  define GFFT_TYPE_C2C  105
#else
#  define GFFT_TYPE_R2C  42
#  define GFFT_TYPE_C2R  44
#  define GFFT_TYPE_C2C  41
#endif
