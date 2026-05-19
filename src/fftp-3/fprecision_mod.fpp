!=================================================================
! MODULES for 3D codes
!
! 2009 Duane Rosenberg and Pablo D. Mininni.
!      National Center for Atmospheric Research.
!=================================================================

!=================================================================

  MODULE fprecision
!
! 
#if defined(GDOUBLE_PRECISION)
      INTEGER, PARAMETER :: GP = KIND(0.0D0)
#elif defined(GSINGLE_PRECISION)
      INTEGER, PARAMETER :: GP = KIND(0.0)
#else
#  error 'MODULE FPRECISION: PRECISION must be GDOUBLE_PRECISION or GSINGLE_PRECISION'
#endif
      SAVE

  END MODULE fprecision
!=================================================================

!=================================================================

  MODULE commtypes
     USE mpi_f08
!
! 
#if defined(GDOUBLE_PRECISION)
      TYPE(mpi_datatype), SAVE :: GC_REAL    = MPI_DOUBLE_PRECISION
      TYPE(mpi_datatype), SAVE :: GC_COMPLEX = MPI_DOUBLE_COMPLEX
#elif defined(GSINGLE_PRECISION)
      TYPE(mpi_datatype), SAVE :: GC_REAL    = MPI_REAL
      TYPE(mpi_datatype), SAVE :: GC_COMPLEX = MPI_COMPLEX
#else
#  error 'MODULE COMMTYPES: PRECISION must be GDOUBLE_PRECISION or GSINGLE_PRECISION'
#endif

  END MODULE commtypes
!=================================================================

