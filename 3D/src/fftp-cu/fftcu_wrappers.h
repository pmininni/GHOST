!=================================================================
! Include file name-mangling of CUFFT API names
!
! 2013 Duane Rosenberg 
!      National Center for Computational Science, ORNL
!=================================================================
#define CONCAT(prefix,name)    prefix/**/name
#if defined(GDOUBLE_PRECISION)
#  define GFLOATBYTESZ     8
#  define GCUFFTEXECR2C    cufftExecD2Z
#  define GCUFFTEXECC2R    cufftExecZ2D
#  define GCUFFTEXECC2C    cufftExecZ2Z
#  define GCUFFTDEFR2C     CUFFT_D2Z
#  define GCUFFTDEFC2R     CUFFT_Z2D
#  define GCUFFTDEFC2C     CUFFT_Z2Z
#elif defined(GSINGLE_PRECISION) 
#  define GFLOATBYTESZ     4
#  define GCUFFTEXECR2C    cufftExecR2C
#  define GCUFFTEXECC2R    cufftExecC2R
#  define GCUFFTEXECC2C    cufftExecC2C
#  define GCUFFTDEFR2C     CUFFT_R2C
#  define GCUFFTDEFC2R     CUFFT_C2R
#  define GCUFFTDEFC2C     CUFFT_C2C
#endif
