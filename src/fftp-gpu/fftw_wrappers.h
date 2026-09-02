!=================================================================
! Include file with the precision-dependent names of the FFTW API
!
! 2009 Duane Rosenberg and Pablo D. Mininni.
!      National Center for Atmospheric Research.
!
! The names are spelled out for each precision instead of being
! built with token pasting, which the preprocessors of gfortran,
! nvfortran and flang do not all support.
!=================================================================
#if defined(GDOUBLE_PRECISION)
#define GFFTW_INIT_THREADS       dfftw_init_threads
#define GFFTW_PLAN_WITH_NTHREADS dfftw_plan_with_nthreads
#define GFFTW_PLAN_MANY_DFT_R2C  dfftw_plan_many_dft_r2c
#define GFFTW_PLAN_MANY_DFT_C2R  dfftw_plan_many_dft_c2r
#define GFFTW_PLAN_MANY_DFT      dfftw_plan_many_dft
#define GFFTW_DESTROY_PLAN       dfftw_destroy_plan
#define GFFTW_EXECUTE_DFT_R2C    dfftw_execute_dft_r2c
#define GFFTW_EXECUTE_DFT_C2R    dfftw_execute_dft_c2r
#define GFFTW_EXECUTE_DFT        dfftw_execute_dft
#elif defined(GSINGLE_PRECISION)
#define GFFTW_INIT_THREADS       sfftw_init_threads
#define GFFTW_PLAN_WITH_NTHREADS sfftw_plan_with_nthreads
#define GFFTW_PLAN_MANY_DFT_R2C  sfftw_plan_many_dft_r2c
#define GFFTW_PLAN_MANY_DFT_C2R  sfftw_plan_many_dft_c2r
#define GFFTW_PLAN_MANY_DFT      sfftw_plan_many_dft
#define GFFTW_DESTROY_PLAN       sfftw_destroy_plan
#define GFFTW_EXECUTE_DFT_R2C    sfftw_execute_dft_r2c
#define GFFTW_EXECUTE_DFT_C2R    sfftw_execute_dft_c2r
#define GFFTW_EXECUTE_DFT        sfftw_execute_dft
#endif
