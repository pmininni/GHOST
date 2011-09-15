!=================================================================
! MODULES for FFTP v3
! Parallel Fast Fourier Transform in 2D using CUDA
!
! 2011 Duane L. Rosenberg and Pablo D. Mininni.
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
      USE iso_c_binding
 
      INTEGER, PARAMETER  :: ikind = IKIND_
      INTEGER, PARAMETER  :: csize = CSIZE_
      INTEGER, PARAMETER  :: nstrip = NSTRIP_
      INTEGER, PARAMETER  :: FFTCU_REAL_TO_COMPLEX = -1
      INTEGER, PARAMETER  :: FFTCU_COMPLEX_TO_REAL =  1
      INTEGER, PARAMETER  :: FFTW_REAL_TO_COMPLEX = -1
      INTEGER, PARAMETER  :: FFTW_COMPLEX_TO_REAL =  1
      INTEGER, PARAMETER  :: FFTW_MEASURE =  0
      INTEGER, PARAMETER  :: FFTW_PATIENT =  0
      INTEGER, PARAMETER  :: FFTW_ESTIMATE=  0
      DOUBLE PRECISION    :: ffttime = 0.0
      DOUBLE PRECISION    :: memtime = 0.0
      DOUBLE PRECISION    :: tratime = 0.0
      TYPE FFTPLAN
         COMPLEX(KIND=GP), POINTER, DIMENSION (:,:)    :: ccarr
         COMPLEX(KIND=GP), POINTER, DIMENSION (:,:)    :: carr
         REAL   (KIND=GP), POINTER, DIMENSION (:,:)    :: rarr
         TYPE     (C_PTR)                              :: cu_ccd_, cu_cd_, cu_rd_
         TYPE     (C_PTR)                              :: pccarr_, pcarr_, prarr_
         INTEGER  (C_INT)                              :: icuplanr_, icuplanc_
         INTEGER :: n
         INTEGER :: szccd_, szcd_, szrd_
         INTEGER, DIMENSION (:), POINTER :: itype1, itype2
      END TYPE FFTPLAN
      SAVE

  END MODULE fftplans
!=================================================================

  MODULE mpivars
!     INCLUDE 'mpif.h'
      INTEGER, SAVE :: ista,iend
      INTEGER, SAVE :: jsta,jend
      INTEGER, SAVE :: ksta,kend
      INTEGER, SAVE :: nprocs,myrank
      INTEGER, SAVE :: ierr

  END MODULE mpivars

!=================================================================
  MODULE cutypes
      IMPLICIT NONE

      
      INTEGER, PARAMETER :: CUFFT_R2C=X'2a'
      INTEGER, PARAMETER :: CUFFT_C2R=X'2c'
      INTEGER, PARAMETER :: CUFFT_C2C=X'29'
      INTEGER, PARAMETER :: CUFFT_D2Z=X'6a'
      INTEGER, PARAMETER :: CUFFT_Z2D=X'6c'
      INTEGER, PARAMETER :: CUFFT_Z2Z=X'69'
      INTEGER, PARAMETER :: cudaHostAllocDefault =0
      INTEGER, PARAMETER :: cudaHostAllocPortable=1
      INTEGER, PARAMETER :: cudaHostAllocMapped  =2
      INTEGER, PARAMETER :: cudaDeviceMapHost    =3
     

  END MODULE cutypes
!=================================================================
