!=================================================================
! Include file name-mangling of FFTW API names
!
! 2009 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!=================================================================
#define CONCAT(prefix,name)    prefix/**/name
#if defined(GDOUBLE_PRECISION)
#define GPMANGLE(name)    CONCAT(dfftw_,name)
#elif defined(GSINGLE_PRECISION) 
#define GPMANGLE(name)    CONCAT(sfftw_,name)
#endif

