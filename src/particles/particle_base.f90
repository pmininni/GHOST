! ===================================================================
! NAME       : particle_base.f90
! DESCRIPTION: Forms base class for all particles
!
! DATE       : 19/02/26 (PDM)
! ===================================================================

module particlebase_mod
  use class_GWorkspace3D
  use gpstate_mod
  use iovar

  ! ================= Base class for all particles ==================
  ! Define an abstract base class
  type, abstract :: ParticleBase
      type(GWorkspace), pointer           :: workspace_
      integer                             :: myrank_   ! MPI rank
      integer                             :: nprocs_   ! MPI procs
      integer                             :: POSITION  ! start of position sector
      integer                             :: nd_       ! problem dimension
      integer                             :: nc_       ! # vector field components
      character(len=8), allocatable       :: sstate_(:)! state member names
      character(len=128)                  :: infile_   ! config file name
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus_
      INTEGER                             :: inittype_
      INTEGER                             :: iinterp_
      INTEGER                             :: iexchtype_
      INTEGER                             :: iouttype_
      INTEGER                             :: intacc_
      INTEGER                             :: wrtunit_
      INTEGER                             :: bcollective_
      INTEGER                             :: itimetype_
      TYPE(GPartComm)                     :: gpcomm_
      TYPE(GPSplineInt)                   :: intop_
      INTEGER                             :: intorder_,itorder_,nd_(3)
      INTEGER                             :: libnds_(3,2),tibnds_(3,2)
      INTEGER                             :: htimers_(GPMAXTIMERS)
      INTEGER                             :: ierr_,iseed_,istep_
      INTEGER                             :: maxparts_,nparts_,npartsm_,nvdb_
      INTEGER                             :: partbuff_,partchunksize_
      INTEGER                             :: comm_
      INTEGER      , ALLOCATABLE, DIMENSION  (:) :: id_,idm_,tmpint_
      INTEGER      , ALLOCATABLE, DIMENSION  (:) :: fpid_
      REAL(KIND=GP), ALLOCATABLE, DIMENSION  (:) :: px_ ,py_ ,pz_
      REAL(KIND=GP), ALLOCATABLE, DIMENSION  (:) :: lvx_,lvy_,lvz_
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:) :: ptmp0_,ptmp1_,ptmp2_,vdb_
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:) :: vk0_,vk1_,vk2_,xk1_
      REAL(KIND=GP)                              :: lxbnds_(3,2),gext_(3)
      REAL(KIND=GP)                              :: delta_(3),invdel_(3)
      CHARACTER(len=1024)                        :: seedfile_,sfile_
      CHARACTER(len=MPI_MAX_ERROR_STRING)        :: serr_
    contains
      procedure(GPart_ctor_interface), deferred :: GPart_ctor => null()
      procedure(Init_interface)      , deferred :: init_particles => null()
      procedure(dvdt_interface)      , deferred :: dvdt_particles => null()
      procedure(write_interface),      deferred :: write_particles => null()
      procedure(state_size_interface), deferred :: state_size  ! Number of states
  END TYPE ParticleBase

  type, abstract, extends(ParticleBase)      :: VelocParticleBase
      integer :: VELOCITY    ! start of velocity sector
  end type VelocParticleBase
   
  type, abstract, extends(VelocParticleBase) :: ChargedParticleBase
  end type VelocParticleBase

  abstract interface
  end interface

CONTAINS

  ! ===================================================================
  ! Concrete methods inherited by all solvers
  ! ===================================================================
  
end module particlebase_mod
