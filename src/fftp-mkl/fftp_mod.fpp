!=================================================================
! MODULES for FFTP v3
! Parallel Fast Fourier Transform in 2D and 3D
!
! 2007 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.uba.ar 
!=================================================================

!=================================================================

  MODULE fftplans
!
! Set the variable ikind to:  4 in 32 bits machines
!                             8 in 64 bits machines
! Set the variable csize to:  8 if L1 cache is <= 64 kb
!                            16 if L1 cache is 128 kb
! The variable nstrip controls strip mining during the 
! transposition. Often set to 1.
!
      USE fprecision
 
      INTEGER, PARAMETER  :: ikind = IKIND_
      INTEGER, PARAMETER  :: csize = CSIZE_
      INTEGER, PARAMETER  :: nstrip = NSTRIP_
      INTEGER             :: dfterr
      INTEGER, PARAMETER  :: FFTW_FORWARD = -1
      INTEGER, PARAMETER  :: FFTW_BACKWARD = +1
      INTEGER, PARAMETER  :: FFTW_REAL_TO_COMPLEX = FFTW_FORWARD
      INTEGER, PARAMETER  :: FFTW_COMPLEX_TO_REAL = FFTW_BACKWARD
      INTEGER, PARAMETER  :: FFTW_MEASURE = 0
      INTEGER, PARAMETER  :: FFTW_EXHAUSTIVE = 8
      INTEGER, PARAMETER  :: FFTW_PATIENT = 32
      INTEGER, PARAMETER  :: FFTW_ESTIMATE = 64
      INTEGER, PARAMETER  :: FFTW_ESTIMATE_PATIENT = 128
      INTEGER             :: hcom,hfft,hmem,htra,htot
      DOUBLE PRECISION    :: comtime = 0.0
      DOUBLE PRECISION    :: ffttime = 0.0
      DOUBLE PRECISION    :: tratime = 0.0
      DOUBLE PRECISION    :: tottime = 0.0
      TYPE FFTPLAN
         TYPE(DFTI_DESCRIPTOR), POINTER :: planr,planc
         INTEGER :: nx,ny,nz
         INTEGER, DIMENSION (:), POINTER :: itype1, itype2
      END TYPE FFTPLAN
      SAVE

  END MODULE fftplans
!=================================================================

  MODULE threads
!
! nth: number of threads used in OpenMP loops and FFTs
!$    USE omp_lib
      INTEGER                            :: nth
      LOGICAL                            :: fft_threads = .false.
      SAVE

  END MODULE threads
!=================================================================


  MODULE mpivars
      INTEGER, SAVE :: ista,iend
      INTEGER, SAVE :: jsta,jend
      INTEGER, SAVE :: ksta,kend
      INTEGER, SAVE :: nprocs,myrank
      INTEGER, SAVE :: provided
      INTEGER, SAVE :: ierr

  END MODULE mpivars
!=================================================================
