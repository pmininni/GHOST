!=================================================================
! Include file name-mangling of CUFFT API names
!
! 2013 Duane Rosenberg 
!      National Center for Computational Science, ORNL
!=================================================================
#define CONCAT(prefix,name)    prefix/**/name
#if defined(GDOUBLE_PRECISION)
#  define GFLOATBYTESZ      8
#  define GMANGLER2C    cufftExecD2Z
#  define GMANGLEC2R    cufftExecZ2D
#  define GMANGLEC2C    cufftExecZ2Z
#  define GDEFR2C       CUFFT_D2Z
#  define GDEFC2R       CUFFT_Z2D
#  define GDEFC2C       CUFFT_Z2Z
#elif defined(GSINGLE_PRECISION) 
#  define GFLOATBYTESZ      4
#  define GMANGLER2C    cufftExecR2C
#  define GMANGLEC2R    cufftExecC2R
#  define GMANGLEC2C    cufftExecC2C
#  define GDEFR2C       CUFFT_R2C
#  define GDEFC2R       CUFFT_C2R
#  define GDEFC2C       CUFFT_C2C
#endif
