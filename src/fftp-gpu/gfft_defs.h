!=================================================================
! Include file with the names of the device FFT library (cuFFT or
! hipFFT) for the selected vendor and precision. The two libraries
! share the C API and the transform type codes; only the prefix and
! the handle type differ (an int for cuFFT, a pointer for hipFFT).
! Names are spelled out instead of built by token pasting, which the
! preprocessors of gfortran, nvfortran and flang do not all support.
!
! 2026 Pablo D. Mininni
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!=================================================================
#if defined(GPU_NVIDIA)
#  define GFFT_HANDLE_T        INTEGER(C_INT)
#  define C_GFFT_PLAN_MANY     "cufftPlanMany"
#  define C_GFFT_DESTROY       "cufftDestroy"
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
#  if defined(GDOUBLE_PRECISION)
#    define C_GFFT_EXEC_R2C    "hipfftExecD2Z"
#    define C_GFFT_EXEC_C2R    "hipfftExecZ2D"
#    define C_GFFT_EXEC_C2C    "hipfftExecZ2Z"
#  else
#    define C_GFFT_EXEC_R2C    "hipfftExecR2C"
#    define C_GFFT_EXEC_C2R    "hipfftExecC2R"
#    define C_GFFT_EXEC_C2C    "hipfftExecC2C"
#  endif
#else
#  error 'FFTP-GPU: GPU_NVIDIA or GPU_AMD must be defined (build with P_GPU)'
#endif

! Transform type codes (identical in cuFFT and hipFFT) and directions
#if defined(GDOUBLE_PRECISION)
#  define GFFT_TYPE_R2C  106
#  define GFFT_TYPE_C2R  108
#  define GFFT_TYPE_C2C  105
#else
#  define GFFT_TYPE_R2C  42
#  define GFFT_TYPE_C2R  44
#  define GFFT_TYPE_C2C  41
#endif
#define GFFT_FORWARD  -1
#define GFFT_INVERSE   1
