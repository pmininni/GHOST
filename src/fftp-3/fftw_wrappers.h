!=================================================================
! Include file name-mangling of FFTW API names
!
! 2009 Duane Rosenberg and Pablo D. Mininni.
!      National Center for Atmospheric Research.
!=================================================================
#define CONCAT(prefix,name)    prefix/**/name
#if defined(GDOUBLE_PRECISION)
#define GPMANGLE(name)    CONCAT(dfftw_,name)
#elif defined(GSINGLE_PRECISION) 
#define GPMANGLE(name)    CONCAT(sfftw_,name)
#endif

