!=================================================================
! Include file name-mangling of CUFFT API names
!
! 2013 Duane Rosenberg 
!      National Center for Computational Science, ORNL
!=================================================================
#define CONCAT(prefix,name)    prefix/**/name
#define cudaErrChk() cfileerr=__FILE__; call scudaErrChk(C_LOC(cfileerr), __LINE__)


#if defined(GDOUBLE_PRECISION)
#  define GFLOATBYTESZ     8
#  define GCUFFTEXECOFFR2C cufftExecOffD2Z
#  define GCUFFTEXECOFFC2R cufftExecOffZ2D
#  define GCUFFTEXECOFFC2C cufftExecOffZ2Z
#  define GCUFFTEXECR2C    cufftExecD2Z
#  define GCUFFTEXECC2R    cufftExecZ2D
#  define GCUFFTEXECC2C    cufftExecZ2Z
#  define GCUFFTDEFR2C     CUFFT_D2Z
#  define GCUFFTDEFC2R     CUFFT_Z2D
#  define GCUFFTDEFC2C     CUFFT_Z2Z
#elif defined(GSINGLE_PRECISION) 
#  define GFLOATBYTESZ     4
#  define GCUFFTEXECOFFR2C cufftExecOffR2C
#  define GCUFFTEXECOFFC2R cufftExecOffC2R
#  define GCUFFTEXECOFFC2C cufftExecOffC2C
#  define GCUFFTEXECR2C    cufftExecR2C
#  define GCUFFTEXECC2R    cufftExecC2R
#  define GCUFFTEXECC2C    cufftExecC2C
#  define GCUFFTDEFR2C     CUFFT_R2C
#  define GCUFFTDEFC2R     CUFFT_C2R
#  define GCUFFTDEFC2C     CUFFT_C2C
#endif
