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
      INCLUDE 'fftw3.f'
 
      INTEGER, PARAMETER  :: ikind = IKIND_
      INTEGER, PARAMETER  :: csize = CSIZE_
      INTEGER, PARAMETER  :: nstrip = NSTRIP_
      INTEGER, PARAMETER  :: FFTW_REAL_TO_COMPLEX = FFTW_FORWARD
      INTEGER, PARAMETER  :: FFTW_COMPLEX_TO_REAL = FFTW_BACKWARD
      INTEGER             :: hcom,hfft,hmem,htra,htot
      DOUBLE PRECISION    :: comtime = 0.0
      DOUBLE PRECISION    :: ffttime = 0.0
      DOUBLE PRECISION    :: tratime = 0.0
      DOUBLE PRECISION    :: tottime = 0.0
      TYPE FFTPLAN
         COMPLEX(KIND=GP), DIMENSION (:,:,:), POINTER :: ccarr
         COMPLEX(KIND=GP), DIMENSION (:,:,:), POINTER :: carr
         REAL(KIND=GP), DIMENSION (:,:,:), POINTER    :: rarr
         INTEGER(kind=ikind) :: planr,planc
         INTEGER :: nx,ny,nz
         INTEGER :: ista,iend,ksta,kend
         INTEGER :: comm,myrank,nprocs
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
!$      IF (ndev.eq.0) THEN
!$         PRINT *, "No devices found in the computing node"
!$      ELSE
!$        tdev = MODULO(myrank,ndev)
!$        IF (ndev.ne.hdev) THEN
!$          PRINT *, "The IDs of the host and target devices do not  &
!$                    conform with the GHOST standard. Please check  &
!$                    that the computing nodes have accelerators, or &
!$                    check whether the numbering of devices in this &
!$                    system goes from 0 to num_devices-1, with the  &
!$                    host device ID equal to the number of devices."
!$          STOP
!$        ENDIF
!$      ENDIF
#endif
      END SUBROUTINE init_offload
    
  END MODULE offloading
!=================================================================

  MODULE mpivars
      INTEGER, SAVE :: ista,iend
      INTEGER, SAVE :: jsta,jend
      INTEGER, SAVE :: ksta,kend
      INTEGER, SAVE :: itsta,itend ! truncated versions
      INTEGER, SAVE :: jtsta,jtend
      INTEGER, SAVE :: ktsta,ktend
      INTEGER, SAVE :: nprocs,myrank
      INTEGER, SAVE :: ntprocs,mytrank
      INTEGER, SAVE :: provided
      INTEGER, SAVE :: ierr

  END MODULE mpivars
!=================================================================
