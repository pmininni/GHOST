!=================================================================
! MODULES for FFTP v2
! Parallel Fast Fourier Transform in 2D and 3D
!
! 2003 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
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
      INCLUDE 'fftw_f77.i'
      INTEGER, PARAMETER  :: ikind = IKIND_
      INTEGER, PARAMETER  :: csize = CSIZE_
      INTEGER, PARAMETER  :: nstrip = NSTRIP_
      DOUBLE PRECISION    :: comtime = 0.0D0
      DOUBLE PRECISION    :: ffttime = 0.0D0
      DOUBLE PRECISION    :: tratime = 0.0D0
      TYPE FFTPLAN
         INTEGER(kind=ikind) :: planr,planc
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
      INTEGER :: nth
      SAVE

  END MODULE threads
!=================================================================

  MODULE offloading
!
! Module with variables for offloading to GPUs using OpenMP. GHOST
! assumes the number of MPI jobs in each node is equal to the
! number of GPUs available in the node. The user must ensure this
! condition is fulfilled.
!$    USE threads
      INTEGER                            :: numdev
      INTEGER                            :: hostdev,targetdev
      LOGICAL                            :: offload
      SAVE
    CONTAINS
      SUBROUTINE init_offload(myrank,ndev,hdev,tdev)
      ! Binds each MPI job to only one GPU target device
        INTEGER, INTENT(IN)    :: myrank
        INTEGER, INTENT(INOUT) :: ndev,hdev,tdev
        ndev = 0
        hdev = 0
        tdev = 0
#if defined(DO_HYBRIDoffl)
!$      ndev = omp_get_num_devices()
!$      hdev = omp_get_initial_device()
!$      tdev = MODULO(myrank,ndev) + 1
#endif
  END MODULE offloading
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
