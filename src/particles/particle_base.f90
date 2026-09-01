! ===================================================================
! NAME       : particle_base.f90
! DESCRIPTION: Forms base class for all particles
!
! DATE       : 19/02/26 (PDM)
! ===================================================================

module particlebase_mod
  use class_GWorkspace3D
  use class_GPSplineInt
  use class_GPartComm
  use gstate_mod
  use commtypes
  use gtimer
  implicit none

  ! ================= Global parameters =============================
  INTEGER,PARAMETER,PUBLIC    :: GPINTRP_CSPLINE=0
  INTEGER,PARAMETER,PUBLIC    :: GPINTRP_LAGINT =1
  INTEGER,PARAMETER,PUBLIC    :: GPEXCHTYPE_NN  =0
  INTEGER,PARAMETER,PUBLIC    :: GPEXCHTYPE_VDB =1
  INTEGER,PARAMETER,PUBLIC    :: GPTIME_STEP    =1
  INTEGER,PARAMETER,PUBLIC    :: GPTIME_COMM    =2
  INTEGER,PARAMETER,PUBLIC    :: GPTIME_SPLINE  =3
  INTEGER,PARAMETER,PUBLIC    :: GPTIME_TRANSP  =4
  INTEGER,PARAMETER,PUBLIC    :: GPTIME_DATAEX  =5
  INTEGER,PARAMETER,PUBLIC    :: GPTIME_INTERP  =6
  INTEGER,PARAMETER,PUBLIC    :: GPTIME_PUPDATE =7
  INTEGER,PARAMETER,PUBLIC    :: GPTIME_GPREAD  =8
  INTEGER,PARAMETER,PUBLIC    :: GPTIME_GPWRITE =9
  INTEGER,PARAMETER,PUBLIC    :: GPSWIPERATE    =100
  INTEGER,PARAMETER,PUBLIC    :: GPMAXTIMERS    =9  ! no. GPTIME parameters
  CHARACTER(len=8),PUBLIC     :: lgext              ! string to hold time index
  CHARACTER(len=6),PUBLIC     :: lgfmtext='(i8.8)'  ! file time index format

  ! ================= Base class for all particles ==================
  ! Define an abstract base class
  type, abstract :: ParticleBase
      type(GWorkspace), pointer           :: workspace_ => null()
      integer                             :: myrank_     ! MPI rank
      integer                             :: nprocs_     ! MPI procs
      integer                             :: POSITION    ! start of position sector
      integer                             :: ndim_       ! problem dimension
      integer                             :: nc_         ! # vector field components
      character(len=8)                    :: sstate_pos_ ! state name of positions
      character(len=8)                    :: sstate_lag_ ! state name of lag velocities
      character(len=128)                  :: infile_     ! config file name
      character(len=128)                  :: odir_,idir_ ! internal class I/O dirs
      INTEGER                             :: iinterp_
      INTEGER                             :: iexchtype_
      INTEGER                             :: iouttype_
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
      INTEGER                             :: partbuff_,partchunksize_,stepcounter_
      TYPE(MPI_Comm)                      :: comm_
      TYPE(MPI_Status)                    :: istatus_
      INTEGER      , ALLOCATABLE, DIMENSION  (:) :: id_,idm_,tmpint_
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:) :: vdb_,ptmp0_,gptmp0_
      REAL(KIND=GP), pointer    , DIMENSION  (:) :: px_  => null()
      REAL(KIND=GP), pointer    , DIMENSION  (:) :: py_  => null()
      REAL(KIND=GP), pointer    , DIMENSION  (:) :: pz_  => null()
      REAL(KIND=GP), pointer    , DIMENSION  (:) :: lvx_ => null()
      REAL(KIND=GP), pointer    , DIMENSION  (:) :: lvy_ => null()
      REAL(KIND=GP), pointer    , DIMENSION  (:) :: lvz_ => null()
      REAL(KIND=GP)                       :: lxbnds_(3,2),gext_(3)
      REAL(KIND=GP)                       :: delta_(3),invdel_(3)
      CHARACTER(len=1024)                 :: seedfile_,sfile_
      CHARACTER(len=MPI_MAX_ERROR_STRING) :: serr_
      logical                             :: hasfeedback_ ! Has feedback on fluid
    contains
      procedure(part_ctor_interface) , deferred :: part_ctor
      procedure(init_interface)      , deferred :: init
      procedure(dpdt_interface)      , deferred :: dpdt
      procedure(end_stage_interface) , deferred :: end_stage
      procedure(feedback_interface)  , deferred :: feedback
      procedure(write_interface),      deferred :: write_pstate
      procedure(state_size_interface), deferred :: state_size  ! Number of states
      procedure, public                         :: SetRandSeed
      procedure, public                         :: GetRandSeed
      procedure, public                         :: AssignLagPos
      procedure, public                         :: io_write_pdb
      procedure, public                         :: io_write_vec
      procedure, public                         :: io_write_euler
      procedure, public                         :: binary_write_lag_co
      procedure, public                         :: binary_write_lag_t0
      procedure, public                         :: ascii_write_lag
      procedure, public                         :: io_read
      procedure, public                         :: binary_read_id_co
      procedure, public                         :: binary_read_pdb_co
      procedure, public                         :: binary_read_pdb_t0
      procedure, public                         :: ascii_read_pdb
      procedure, public                         :: io_readvec
      procedure, public                         :: io_readvec_scatter_
      procedure, public                         :: EulerToLag
      procedure, public                         :: MakePeriodicP
      procedure, public                         :: MakePeriodicZ
      procedure, public                         :: MakePeriodicExt
      procedure, public                         :: Part_Delete
      procedure, public                         :: GetLocalWrk
      procedure, public                         :: GetLocalWrk_aux
      procedure, public                         :: CopyLocalWrk
      procedure, public                         :: GetVDB
      procedure, public                         :: GetVel
      procedure, public                         :: GetTime
      procedure, public                         :: GetNParts
      procedure, public                         :: GetLoadBal
      procedure, public                         :: PartNumConsistent
      procedure, public                         :: R3toR3
      procedure, public                         :: ResizeArrays
  end type ParticleBase

  type, abstract, extends(ParticleBase)      :: VelocParticleBase
      integer                             :: VELOCITY    ! start of velocity sector
      character(len=8)                    :: sstate_vel_ ! state name of velocity
  end type VelocParticleBase

  abstract interface
     subroutine part_ctor_interface(this, infile, pde, workspace, pstate, pstate_cpy)
       use class_GWorkspace3D
       use equationbase_mod
       use gpstate_mod
       import :: ParticleBase
       class(ParticleBase), intent(inout)              :: this
       class(EquationBase), intent   (in)              :: pde
       type   (GWorkspace), intent(inout), target      :: workspace
       type  (GPStateComp), intent(inout), allocatable :: pstate(:), pstate_cpy(:)
       character   (len=*), intent   (in)              :: infile
     end subroutine part_ctor_interface

     subroutine init_interface(this)
       import :: ParticleBase
       class(ParticleBase),         intent(inout) :: this       
     end subroutine init_interface

     subroutine dpdt_interface(this, time, pde, fluidstate, pstate, dt, dpdtout)
       use equationbase_mod
       use gpstate_mod
       import :: ParticleBase
       class(ParticleBase),         intent(inout) :: this
       class(EquationBase),         intent   (in) :: pde
       real      (kind=GP),         intent   (in) :: time, dt
       type   (GStateComp),         intent   (in) :: fluidstate(:)
       type  (GPStateComp), target, intent   (in) :: pstate(:) 
       type  (GPStateComp),         intent(inout) :: dpdtout(:) 
     end subroutine dpdt_interface

     subroutine end_stage_interface(this, upin, upout)
       use gpstate_mod
       import :: ParticleBase
       class(ParticleBase), intent(inout)         :: this
       type  (GPStateComp), intent(inout)         :: upin (:)
       type  (GPStateComp), intent(inout)         :: upout(:)
     end subroutine end_stage_interface

     subroutine feedback_interface(this, pstate, feedback)
       use gstate_mod
       use gpstate_mod
       import :: ParticleBase
       class(ParticleBase), intent   (in)         :: this
       type  (GPStateComp), intent   (in)         :: pstate(:)
       type   (GStateComp), intent(inout)         :: feedback(:)
     end subroutine feedback_interface

     subroutine write_interface(this, time, pde, fluidstate, pstate)
       use equationbase_mod
       use gpstate_mod
       import :: ParticleBase
       class(ParticleBase),         intent(inout) :: this
       class(EquationBase),         intent   (in) :: pde
       real      (kind=GP),         intent   (in) :: time
       type   (GStateComp),         intent   (in) :: fluidstate(:)
       type  (GPStateComp), target, intent   (in) :: pstate(:) 
     end subroutine write_interface

     function state_size_interface(this) result(num)
       import :: ParticleBase
       class(ParticleBase),         intent   (in) :: this
       integer                                    :: num
     end function state_size_interface
  end interface

CONTAINS

  ! ===================================================================
  ! Concrete methods inherited by all solvers
  ! ===================================================================
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : SetRandSeed
  !!  DESCRIPTION: Sets seed for random number generator
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance
  !!    seed    : value of seed
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SetRandSeed(this, iseed)
    USE random
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT) :: this
    INTEGER             ,INTENT   (IN) :: iseed
    this%iseed_     = iseed;
    CALL prandom_seed(this%iseed_) 
  END SUBROUTINE SetRandSeed

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : GetRandSeed
  !!  DESCRIPTION: Gets random seed
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION GetRandSeed(this) result(get_res)
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT) :: this
    INTEGER                            :: get_res
    get_res = this%iseed_
  END FUNCTION GetRandSeed


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : AssignLagPos
  !!  DESCRIPTION: Assigns px,py,pz pointers to a particle state.
  !!               Must be called before any routine that needs
  !!               access to this%px_, this%py_, this%pz_.
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance
  !!    pstate  : A particle state
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE AssignLagPos(this, pstate)
    use gpstate_mod
    CLASS(ParticleBase) ,INTENT(INOUT)      :: this
    type  (GPStateComp), intent(in), target :: pstate(:)
    nullify (this%px_)
    this%px_ => pstate(this%POSITION  )%rcomp
    nullify (this%py_)
    this%py_ => pstate(this%POSITION+1)%rcomp
    nullify (this%pz_)
    this%pz_ => pstate(this%POSITION+2)%rcomp
  END SUBROUTINE AssignLagPos


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : io_write_pdb
  !!  DESCRIPTION: Does write of Lagrangian position d.b. to file. 
  !!               Position of the particle structure in file is the
  !!               particle's id. Main entry point for both ASCII and
  !!               binary writes. This is a special API for outputting 
  !!               internal part. d.b. using underlying mehods that
  !!               handle binary (collective or not) and ASCII I/O.
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance
  !!    iunit   : unit number
  !!    dir     : output directory
  !!    spref   : filename prefix
  !!    nmd     : time index
  !!    time    : real time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE io_write_pdb(this, iunit, dir, spref, nmb, time)
    USE fprecision
    USE commtypes
    USE mpivars
    USE grid
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT) :: this
    REAL(KIND=GP),INTENT   (IN)        :: time
    REAL(KIND=GP)                      :: prec(3)
    INTEGER,INTENT(IN)                 :: iunit
    INTEGER                            :: ht,j,nt
    INTEGER(kind=MPI_OFFSET_KIND)      :: offset
    TYPE(MPI_File)                     :: fh
    CHARACTER(len=*),INTENT(IN)        :: dir
    CHARACTER(len=*),INTENT(IN)        :: nmb
    CHARACTER(len=*),INTENT(IN)        :: spref

    ! Do a sanity check:
    !!  CALL MPI_ALLREDUCE(this%nparts_,nt,1,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
    !!  IF ( nt .NE. this%maxparts_ ) THEN
    !!    WRITE(*,*) this%myrank_, ': io_write_pdb: particle inconsistency: no. required=',&
    !!               this%maxparts_,' no. found=',nt
    !!    STOP
    !!  ENDIF

    IF (this%iexchtype_.EQ.GPEXCHTYPE_NN) THEN
      IF (this%bcollective_.EQ.0) THEN
        CALL this%gpcomm_%VDBSynch_t0(this%ptmp0_,this%maxparts_,this%id_,&
             this%px_,this%py_,this%pz_,this%nparts_)
      END IF
    ELSE IF (this%iexchtype_.EQ.GPEXCHTYPE_VDB) THEN
      CALL this%gpcomm_%VDBSynch(this%vdb_,this%maxparts_,this%id_, &
                     this%px_,this%py_,this%pz_,this%nparts_,this%ptmp0_)
      ! Store global VDB data into temp array:
!$omp parallel do
      DO j = 1, this%maxparts_
        this%ptmp0_(1,j) = this%vdb_(1,j)
        this%ptmp0_(2,j) = this%vdb_(2,j)
        this%ptmp0_(3,j) = this%vdb_(3,j)
      ENDDO
    ENDIF

    CALL GTStart(this%htimers_(GPTIME_GPWRITE))
    CALL GTInitHandle(ht,GT_WTIME)

    IF ( this%iouttype_ .EQ. 0 ) THEN
      IF ( this%bcollective_ .EQ. 1 ) THEN
        ! pass in the current linear _local_ particle coord arrays
        IF ( this%wrtunit_ .EQ. 1 ) THEN ! rescale coordinates to box units
           this%ptmp0_(1,:) = this%px_*this%delta_(1)
           this%ptmp0_(2,:) = this%py_*this%delta_(2)
           this%ptmp0_(3,:) = this%pz_*this%delta_(3)
           CALL this%binary_write_lag_co(iunit,dir,spref,nmb,time,this%nparts_, &
                this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:))
        ELSE
           CALL this%binary_write_lag_co(iunit,dir,spref,nmb,time,this%nparts_, &
                this%px_,this%py_,this%pz_)
        ENDIF
      ELSE
        ! pass in the synched-up VDB (copied to ptmp0_):
        IF ( this%wrtunit_ .EQ. 1 ) THEN ! rescale coordinates to box units
           this%ptmp0_(1,:) = this%ptmp0_(1,:)*this%delta_(1)
           this%ptmp0_(2,:) = this%ptmp0_(2,:)*this%delta_(2)
           this%ptmp0_(3,:) = this%ptmp0_(3,:)*this%delta_(3)
        ENDIF
        CALL this%binary_write_lag_t0(iunit,dir,spref,nmb,time,this%maxparts_,  &
                this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:))
      ENDIF
    ELSE
      ! pass in the synched-up VDB (copied to ptmp0_):
      IF ( this%wrtunit_ .EQ. 1 ) THEN ! rescale coordinates to box units
         this%ptmp0_(1,:) = this%ptmp0_(1,:)*this%delta_(1)
         this%ptmp0_(2,:) = this%ptmp0_(2,:)*this%delta_(2)
         this%ptmp0_(3,:) = this%ptmp0_(3,:)*this%delta_(3)
      ENDIF
      CALL this%ascii_write_lag(iunit,dir,spref,nmb,time,this%maxparts_, &
           this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:))
    ENDIF

    CALL GTAcc(this%htimers_(GPTIME_GPWRITE))
    CALL GTStop(ht)
    if(this%myrank_.eq.0) write(*,*)'io_write_pdb: file: ', spref,'  write time: ', GTGetTime(ht)
    CALL GTFree(ht)
  END SUBROUTINE io_write_pdb

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : io_write_vec
  !!  DESCRIPTION: Does write of Lagrangian vector that is
  !!               currently stored. This vector may not be the
  !!               advecting velocity if a call to SetLagVec
  !!               is made with a different vector. This is a 
  !!               special API for outputting the class' internal 
  !!               'velocity' vector using underlying mehods that
  !!               handle binary (collective or not) and ASCII I/O.
  !!               
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance
  !!    iunit   : unit number
  !!    dir     : output directory
  !!    spref   : filename prefix
  !!    nmd     : time index
  !!    time    : real time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE io_write_vec(this, iunit, dir, spref, nmb, time)
    USE fprecision
    USE commtypes
    USE mpivars
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT) :: this
    REAL(KIND=GP),INTENT   (IN)        :: time
    INTEGER,INTENT(IN)                 :: iunit
    CHARACTER(len=*),INTENT(IN)        :: dir
    CHARACTER(len=*),INTENT(IN)        :: nmb
    CHARACTER(len=*),INTENT(IN)        :: spref
    INTEGER                            :: gsum,ht,j

    IF ( .NOT.this%PartNumConsistent(this%nparts_,gsum) ) THEN
      write(*,*)'io_write_vec: global sum=',gsum,' maxparts=',this%maxparts_
      IF ( this%myrank_.eq.0 ) THEN
        WRITE(*,*) 'io_write_vec: Inconsistent particle count'
        STOP
      ENDIF
    ENDIF

    CALL GTInitHandle(ht,GT_WTIME)
    ! If doing non-collective binary or ascii writes, synch up vector:
    IF ((this%iouttype_.EQ.0 .AND. this%bcollective_.EQ.0).OR.this%iouttype_.EQ.1 ) THEN
      CALL this%gpcomm_%VDBSynch_t0(this%gptmp0_,this%maxparts_,this%id_,       &
                                 this%lvx_,this%lvy_,this%lvz_,this%nparts_)
    ENDIF

    IF ( this%iouttype_ .EQ. 0 ) THEN
       IF ( this%bcollective_.EQ. 1 ) THEN
        CALL this%binary_write_lag_co(iunit,dir,spref,nmb,time, this%nparts_,   &
                                 this%lvx_,this%lvy_,this%lvz_)
      ELSE
        CALL this%binary_write_lag_t0(iunit,dir,spref,nmb,time, this%maxparts_, &
                                 this%gptmp0_(1,:),this%gptmp0_(2,:),this%gptmp0_(3,:));
      ENDIF
    ELSE
      CALL this%ascii_write_lag(iunit,dir,spref,nmb,time,       this%maxparts_, &
                                 this%gptmp0_(1,:),this%gptmp0_(2,:),this%gptmp0_(3,:));
    ENDIF
    CALL GTStop(ht)
    if(this%myrank_.eq.0) write(*,*)'io_write_vec: file: ', spref,'  write time: ', GTGetTime(ht)
    CALL GTFree(ht)
  END SUBROUTINE io_write_vec


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : io_write_euler
  !!  DESCRIPTION: Converts specified Eulerian real-space variable to
  !!               a Lagrangian quantity by interpolating to particle positions;
  !!               does write of Lagrangian variable to file, depending on 
  !!               class settings. This routine overwrites lvx,lvy in the class.
  !!  ARGUMENTS  :
  !!    this     : 'this' class instance
  !!    iunit    : unit number
  !!    dir      : output directory
  !!    fname    : filename prefix
  !!    nmb      : time index
  !!    time     : real time
  !!    evar     : Eulerian data from which to compute Lagrangian 
  !!               quantity: theta(y) = theta(x(y),t). Interpolation
  !!               of evar is done internally before write. Note that
  !!               data in evar is lost on exit.
  !!    doupdate : if true, do interp point update in interpolator; else don't
  !!    tmp1,tmp2: tmp arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE io_write_euler(this, iunit, dir, spref, nmb, time, evar, doupdate, tmp1, tmp2)
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT)           :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:,:,:) :: evar
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:,:,:) :: tmp1,tmp2
    REAL(KIND=GP),INTENT   (IN)                  :: time
    INTEGER      ,INTENT   (IN)                  :: iunit
    INTEGER                                      :: offset,nt,szint,szreal
    INTEGER                                      :: ht,j
    TYPE(MPI_File)                               :: fh
    LOGICAL      ,INTENT   (IN)                  :: doupdate
    CHARACTER(len=*), INTENT(IN)                 :: dir
    CHARACTER(len=*)  , INTENT(IN)               :: nmb
    CHARACTER(len=*)  , INTENT(IN)               :: spref
    CHARACTER(len=1024)                          :: sfile
    logical                                      :: bret

    CALL this%EulerToLag(this%lvy_,this%nparts_,evar,doupdate,tmp1,tmp2)
    CALL GTInitHandle(ht,GT_WTIME)

    ! If doing non-collective binary or ascii writes, synch up vector:
    IF ( this%iouttype_.EQ.0 .AND. this%bcollective_.EQ.0 .OR. this%iouttype_.EQ.1 ) THEN
       CALL this%gpcomm_%LagSynch_t0(this%lvx_,this%maxparts_,this%id_,this%lvy_,this%nparts_)
    ENDIF

    IF ( this%iouttype_ .EQ. 0 ) THEN
      IF ( this%bcollective_.EQ. 1 ) THEN
        CALL this%binary_write_lag_co(iunit,dir,spref,nmb,time, this%nparts_, this%lvy_)
      ELSE
        CALL this%binary_write_lag_t0(iunit,dir,spref,nmb,time,this%maxparts_,this%lvx_)
      ENDIF
    ELSE
      CALL this%ascii_write_lag      (iunit,dir,spref,nmb,time,this%maxparts_,this%lvx_)
    ENDIF
    CALL GTStop(ht)
    if(this%myrank_.eq.0) write(*,*)'io_write_euler: file: ', spref,'  write time: ', GTGetTime(ht)
    CALL GTFree(ht)
  END SUBROUTINE io_write_euler


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : binary_write_lag_co
  !!  DESCRIPTION: Does collective binary write of Lagrangian field to file. 
  !!               Position of the particle structure in file is the
  !!               particle's id. This method allows for up to 3
  !!               Lagranian variables to be outputted. At least one
  !!               variable _must_ be present (fld0). Do not use keywords
  !!               to specify optional arguments. 
  !!               Note that this call will have all MPI tasks write
  !!               their data collectively, so no 'synching' of data 
  !!               is required on input.
  !!
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance
  !!    iunit   : unit number
  !!    dir     : output directory
  !!    spref   : filename prefix
  !!    nmd     : time index
  !!    time    : real time
  !!    fld0-2  : Lagrangian field
  !!    np      : no. particles to write
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE binary_write_lag_co(this, iunit, dir, spref, nmb, time, np, &
             fld0, fld1, fld2)
    USE fprecision
    USE commtypes
    USE mpivars
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT)   :: this
    REAL(KIND=GP),INTENT   (IN)          :: time
    REAL(KIND=GP),INTENT   (IN)          :: fld0(this%maxparts_)
    REAL(KIND=GP),INTENT   (IN),OPTIONAL,DIMENSION(this%maxparts_) :: fld1,fld2
    REAL(KIND=GP)                        :: vout(3)
    INTEGER,INTENT(IN)                   :: iunit
    INTEGER,INTENT(IN)                   :: np
    INTEGER                              :: nerr,nt,nv,szreal
    INTEGER(kind=MPI_OFFSET_KIND)        :: offset
    TYPE(MPI_File)                       :: fh
    CHARACTER(len=128),INTENT(IN)        :: dir
    CHARACTER(len=*),INTENT(IN)          :: nmb
    CHARACTER(len=*),INTENT(IN)          :: spref
    INTEGER                              :: j,gc,lc

    nv = 1 
    IF ( present(fld1) ) nv=nv+1
    IF ( present(fld2) ) nv=nv+1
    ! Must write part. data to correct position in file:
    CALL MPI_TYPE_SIZE(GC_REAL    ,szreal,this%ierr_)
    CALL MPI_FILE_OPEN(this%comm_,trim(dir) // '/' // trim(spref) // &
         '.' // nmb // '.lag',MPI_MODE_CREATE+MPI_MODE_WRONLY, &
          MPI_INFO_NULL,fh,this%ierr_)
    IF ( this%ierr_ .NE. MPI_SUCCESS ) THEN
      CALL MPI_ERROR_STRING(this%ierr_, this%serr_, nerr,ierr);
      WRITE(*,*) 'binary_write_lag_co: Error reading opening : ', trim(dir) // '/' // trim(spref) // &
         '.' // nmb // '.lag: ', trim(this%serr_)
      STOP
    ENDIF
    offset = 0
    CALL MPI_FILE_WRITE_AT_ALL(fh,offset,real(this%maxparts_,kind=GP),1,GC_REAL,this%istatus_,       &
         this%ierr_)
    offset = INT(szreal,MPI_OFFSET_KIND)
    CALL MPI_FILE_WRITE_AT_ALL(fh,offset,time   ,1,GC_REAL,this%istatus_,this%ierr_)
    gc = 0
    DO j = 1, np
      offset  = INT((nv*this%id_(j)+2),MPI_OFFSET_KIND) * INT(szreal,MPI_OFFSET_KIND)
      vout(1) = fld0(j)
      IF ( present(fld1) ) vout(2) = fld1(j)
      IF ( present(fld2) ) vout(3) = fld2(j)
      CALL MPI_FILE_WRITE_AT(fh,offset,vout,nv,GC_REAL,this%istatus_,this%ierr_)
      CALL MPI_GET_COUNT(this%istatus_,GC_REAL,lc,this%ierr_)
      gc = gc+lc
    ENDDO
    CALL MPI_FILE_CLOSE(fh,this%ierr_)
    IF ( gc .NE. np*nv ) THEN
      WRITE(*,*)this%myrank_, &
        ': binary_write_lag_co: insufficient amount of data written; no. required=',&
        np*nv,' no. written=',gc
      STOP
    ENDIF
  END SUBROUTINE binary_write_lag_co


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : binary_write_lag_t0
  !!  DESCRIPTION: Does binary write of Lagrangian field to file writing 
  !!               to task 0. Position of the particle structure in file is the
  !!               particle's id. This method allows for up to 3
  !!               Lagranian variables to be outputted. At least one
  !!               variable _must_ be present (fld0). Do not use keywords
  !!               to specify optional arguments. 
  !!               Note that this call WILL NOT synch up the particle 
  !!               d.b.. The reason is that this is method can be used for
  !!               _any_ Lagrangian field, not just the particles. It is 
  !!               assumed that the Lagrange fields, f0-2 are already 'synched'
  !!               on entry, since only task 0 will write.
  !!
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance
  !!    iunit   : unit number
  !!    dir     : output directory
  !!    spref   : filename prefix
  !!    nmd     : time index
  !!    time    : real time
  !!    fld0-2  : Lagrangian field
  !!    np      : no. particles to write
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE binary_write_lag_t0(this, iunit, dir, spref, nmb, time, np, &
             fld0, fld1, fld2)
    USE fprecision
    USE commtypes
    USE mpivars
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT)   :: this
    REAL(KIND=GP),INTENT   (IN)          :: time
    REAL(KIND=GP),INTENT   (IN)          :: fld0(this%maxparts_)
    REAL(KIND=GP),INTENT   (IN),OPTIONAL,DIMENSION(this%maxparts_) :: fld1,fld2
    INTEGER,INTENT(IN)                   :: iunit
    INTEGER,INTENT(IN)                   :: np
    INTEGER                              :: nerr,nv
    TYPE(MPI_File)                       :: fh
    CHARACTER(len=128),INTENT(IN)        :: dir
    CHARACTER(len=*),INTENT(IN)          :: nmb
    CHARACTER(len=*),INTENT(IN)          :: spref
    INTEGER                              :: j

    IF ( this%myrank_ .EQ. 0 ) THEN
      nv = 1
!$omp parallel do
      DO j = 1, np
        this%ptmp0_(1,j) = fld0(j)
      ENDDO
      IF ( present(fld1) ) THEN
!$omp parallel do
        DO j = 1, np
          this%ptmp0_(2,j) = fld1(j)
        ENDDO
        nv = nv+1
      ENDIF
      IF ( present(fld2) ) THEN
!$omp parallel do
        DO j = 1, np
          this%ptmp0_(3,j) = fld2(j)
        ENDDO
        nv = nv+1
      ENDIF
      ! 'access=stream' is required here:
      OPEN(iunit,file=trim(dir) // '/' // trim(spref) // &
                       '.' // nmb // '.lag',form='unformatted',access='stream',&
                       iostat=this%ierr_)
      IF ( this%ierr_.NE.0 ) THEN
        WRITE(*,*)'binary_write_lag_t0: could not open file for reading: ',&
        trim(dir)// '/' // trim(spref) // '.' // nmb //  '.lag'
        STOP
      ENDIF
      WRITE(iunit) real(np,kind=GP)
      WRITE(iunit) time
      WRITE(iunit) this%ptmp0_(1:nv,1:this%maxparts_)
      CLOSE(iunit)
    ENDIF
  END SUBROUTINE binary_write_lag_t0


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : ascii_write_lag
  !!  DESCRIPTION: Does ASCII write of Lagrangian fld to file.
  !!               The local MPI tasks write to a file with prefix
  !!               spref, in the following format:
  !!                     dir/spref.TTT.PPP.txt
  !!               where TTT is the time index, given by nmb, and
  !!               PPP is the MPI rank.  This method allows for up to 3
  !!               Lagranian variables to be outputted. At least one
  !!               variable _must_ be present (fld0). Do not use keywords
  !!               to specify optional arguments.
  !!               Note that this call will have ony MPI task 0 write
  !!               data, so data should be 'synched' before calling. 
  !!               This is not done here.
  !!
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance
  !!    iunit   : unit number
  !!    dir     : output directory
  !!    spref   : filename prefix
  !!    nmd     : time index
  !!    time    : real time
  !!    fld0-2  : Lagrangian fields
  !!    np      : no. particles to write
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ascii_write_lag(this, iunit, dir, spref, nmb, time, np, &
            fld0, fld1, fld2)
    USE fprecision
    USE commtypes
    USE mpivars
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT) :: this
    REAL(KIND=GP),INTENT   (IN)        :: time
    REAL(KIND=GP),INTENT   (IN)        :: fld0(this%maxparts_)
    REAL(KIND=GP),INTENT   (IN), OPTIONAL, DIMENSION(this%maxparts_) :: fld1,fld2
    REAL(KIND=GP)                      :: vout(3)
    INTEGER,INTENT   (IN)              :: iunit
    INTEGER,INTENT   (IN)              :: np
    INTEGER                            :: j,nv
    CHARACTER(len=*),INTENT(IN)        :: dir
    CHARACTER(len=*),INTENT(IN)        :: nmb
    CHARACTER(len=*),INTENT(IN)        :: spref
    CHARACTER(len=3)                   :: sind

    nv = 1 
    IF ( present(fld1) ) nv=nv+1
    IF ( present(fld2) ) nv=nv+1
    ! Write global VDB, with time header, indexed only
    ! by time index: dir/spref.TTT.txt:
    IF ( this%myrank_.EQ.0 ) THEN
      OPEN(iunit,file=trim(dir)// '/' // trim(spref) // '.' // &
            nmb //  '.txt')
      WRITE(iunit,*) np
      WRITE(iunit,*) time
      DO j = 1, np
        vout(1) = fld0(j)
        IF ( present(fld1) ) vout(2) = fld1(j)
        IF ( present(fld2) ) vout(3) = fld2(j)
        WRITE(iunit,600) vout(1:nv)
  600   FORMAT(3(E23.15,1X))
      ENDDO
     CLOSE(iunit)
   ENDIF
  END SUBROUTINE ascii_write_lag

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : io_read
  !!  DESCRIPTION: Does read of Lagrangian particle data from file,
  !!               and scattering of work to correct MPI tasks. 
  !!               This is the main entry point for both binary and
  !!               ASCII reads.
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance
  !!    pstate  : particle state vector where positions are stored
  !!    iunit   : unit number
  !!    dir     : input directory
  !!    spref   : filename prefix
  !!    nmb     : time index. If == '' (of 0-length), then no name mangling 
  !!              is done, and 'spref' is assumed to be the fully
  !!              resolved filename. 'dir' will be ignored in this case.
  !!    id      : optional. Dummy array for ids
  !!    lx,ly,lz: optional. If specified, loads PDB into these arrays.
  !!              If not specified, loads data into member data arrays.
  !!              All must be specified if one is, as must nl....
  !!    nl      : optional. If specified, copies amount of local work
  !!              (size of lx, ly, lz) to this variable for output. 
  !!              Note that nl must be specified if lx,ly,lz are.
  !!    opiotype: optional. Overrides member data iouttype_ if specified.
  !!    opbcoll : optional. Overrides member data bcollective_ if specified.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE io_read(this,pstate,iunit,dir,spref,nmb,id,lx,ly,lz,nl, &
                     opiotype,opbcoll)
    use gpstate_mod
    USE fprecision
    USE commtypes
    USE mpivars
    USE grid
    IMPLICIT NONE
    CLASS(ParticleBase)    ,INTENT(INOUT)     :: this
    TYPE  (GPStateComp)    ,INTENT(INOUT)     :: pstate(:)
    REAL(KIND=GP),INTENT(OUT),OPTIONAL,   DIMENSION(:) :: lx,ly,lz
    REAL(KIND=GP)                             :: rvar,time
    INTEGER      ,INTENT(OUT),OPTIONAL,   DIMENSION(:) :: id
    INTEGER,INTENT(IN)                        :: iunit
    INTEGER,INTENT(INOUT),OPTIONAL            :: nl
    INTEGER,INTENT(IN),OPTIONAL               :: opbcoll,opiotype
    INTEGER                                   :: j,ng
    INTEGER                                   :: bcoll,iotype
    INTEGER(kind=MPI_OFFSET_KIND)             :: offset
    TYPE(MPI_File)                            :: fh
    CHARACTER(len=128),INTENT   (IN)          :: dir
    CHARACTER(len=*)  ,INTENT   (IN)          :: nmb
    CHARACTER(len=*)  ,INTENT   (IN)          :: spref

    IF ( present(opiotype) ) THEN
      iotype = opiotype
    ELSE
      iotype = this%iouttype_
    ENDIF
    IF ( present(opbcoll) ) THEN
      bcoll = opbcoll
    ELSE
      bcoll = this%bcollective_
    ENDIF

    IF (this%iexchtype_.EQ.GPEXCHTYPE_NN) THEN
      IF (bcoll.EQ.1) THEN
        IF (len_trim(nmb).gt.0 ) THEN
          CALL this%binary_read_id_co(iunit, trim(dir)          &
              // '/' // trim(spref) // '.' // nmb //'.lag')
        ELSE
          CALL this%binary_read_id_co(iunit, trim(spref))
        ENDIF
      ELSE
        IF (len_trim(nmb).gt.0 ) THEN
          CALL this%binary_read_pdb_t0(iunit, trim(dir)         &
              // '/' // trim(spref) // '.' // nmb //'.lag',     &
              time,this%ptmp0_,.true.)
        ELSE
          CALL this%binary_read_pdb_t0(iunit, trim(spref),      &
              time,this%ptmp0_,.true.)
        END IF
      END IF
      IF (size(pstate(1)%rcomp) .NE. this%partbuff_) THEN
        CALL GPState_resize(pstate,this%partbuff_)
      END IF
    END IF

    CALL this%AssignLagPos(pstate)
    
    CALL GTStart(this%htimers_(GPTIME_GPREAD))
    IF ( iotype .EQ. 0 ) THEN   ! Binary files
      IF ( bcoll.EQ. 1 ) THEN   ! collective binary
        IF (len_trim(nmb).gt.0 ) THEN
          CALL this%binary_read_pdb_co(iunit, trim(dir)         &
              // '/' // trim(spref) // '.' // nmb // '.lag',time,this%ptmp0_)
        ELSE
          CALL this%binary_read_pdb_co(iunit, trim(spref),  time,this%ptmp0_)
        ENDIF
      ELSE                      ! master thread binary
        IF (len_trim(nmb).gt.0 ) THEN
          CALL this%binary_read_pdb_t0(iunit, trim(dir)         &
              // '/' // trim(spref) // '.' // nmb // '.lag',time,this%ptmp0_)
        ELSE
          CALL this%binary_read_pdb_t0(iunit, trim(spref),  time,this%ptmp0_)
        ENDIF
      ENDIF
    ELSE                        ! ASCII files
      IF (len_trim(nmb).gt.0 ) THEN
        CALL this%ascii_read_pdb (iunit, trim(dir)              &
              // '/' // trim(spref) // '.' // nmb // '.txt',time,this%ptmp0_)
      ELSE
        CALL this%ascii_read_pdb (iunit,trim(spref),        time,this%ptmp0_)
      ENDIF
    ENDIF
    ! rescale coordinates from box units
    IF (this%wrtunit_ .EQ. 1) THEN
       this%ptmp0_(1,:) = this%ptmp0_(1,:)*this%invdel_(1)
       this%ptmp0_(2,:) = this%ptmp0_(2,:)*this%invdel_(2)
       this%ptmp0_(3,:) = this%ptmp0_(3,:)*this%invdel_(3)
    ENDIF
    CALL GTAcc(this%htimers_(GPTIME_GPREAD))

    IF ( (this%iexchtype_.EQ.GPEXCHTYPE_VDB).AND.               &
    .NOT.(present(id).and.present(lx).and.present(ly).and.present(lz).and.present(nl)) ) THEN 
      ! Store in member data arrays
      CALL this%GetLocalWrk(this%id_,this%px_,this%py_,this%pz_,&
                             this%nparts_,this%ptmp0_,this%maxparts_)
    ELSE IF(this%iexchtype_.EQ.GPEXCHTYPE_VDB) THEN
      ! Store in specified input arrays
      CALL this%GetLocalWrk(id,lx,ly,lz,                        &
                             nl,this%ptmp0_,this%maxparts_)
    ELSE IF (this%iexchtype_.EQ.GPEXCHTYPE_NN) THEN
      DO j = 1,this%nparts_
        this%px_(j) = this%ptmp0_(1,j)
        this%py_(j) = this%ptmp0_(2,j)
        this%pz_(j) = this%ptmp0_(3,j)
      END DO
    END IF

    CALL MPI_ALLREDUCE(this%nparts_,ng,1,MPI_INTEGER,           &
                       MPI_SUM,this%comm_,this%ierr_)
    IF ( this%myrank_.EQ.0 .AND. ng.NE.this%maxparts_ ) THEN
      WRITE(*,*)'io_read: inconsistent d.b.: expected: ',       &
                 this%maxparts_, '; found: ',ng
      STOP
    ENDIF

    IF ( .NOT.(present(lx).and.present(ly).and.present(lz)) ) THEN 
      ! If there is a global VDB for data 'exchanges', create it here,
      ! but only if data is loaded into member arrays:
      IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN
!$omp parallel do
        DO j = 1, this%maxparts_
          this%vdb_(1:3,j) = this%ptmp0_(1:3,j)
        ENDDO
      ENDIF
    ENDIF
  END SUBROUTINE io_read

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : binary_read_id_co
  !!  DESCRIPTION: Does read of binary Lagrangian particle data from file, 
  !!               collectively to determine corresponding ids.
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance
  !!    iunit   : unit number
  !!    sfile   : fully resolved file name 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE binary_read_id_co(this,iunit,sfile)
    USE fprecision
    USE commtypes
    USE mpivars
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT)        :: this
    REAL(KIND=GP)                             :: rvar,time
    INTEGER,INTENT(IN)                        :: iunit
    INTEGER                                   :: i,j,nerr,szreal,nr,nb
    INTEGER(kind=MPI_OFFSET_KIND)             :: offset
    TYPE(MPI_File)                            :: fh
    CHARACTER(len=*),INTENT   (IN)            :: sfile

    CALL MPI_TYPE_SIZE(GC_REAL,szreal,this%ierr_)
    CALL MPI_FILE_OPEN(this%comm_,trim(sfile),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,this%ierr_)
    IF ( this%ierr_ .NE. MPI_SUCCESS ) THEN
      CALL MPI_ERROR_STRING(this%ierr_, this%serr_, nerr,ierr);
      WRITE(*,*) 'binary_read_pdb_count: Error reading opening : ', trim(sfile),& 
                trim(this%serr_)
      STOP
    ENDIF
  
    ! Must read part. data from correct spot in file:
    offset = 0
    CALL MPI_FILE_READ_AT_ALL(fh,offset,rvar,1,GC_REAL,this%istatus_,this%ierr_)    !  no.parts
    IF ( int(rvar).NE.this%maxparts_ ) THEN
      WRITE(*,*) 'binary_read_pdb_count: Attempt to read incorrect number of particles: required:', &
                  this%maxparts_,' no. read: ',int(rvar)
      WRITE(*,*) 'binary_read_pdb_count: Error reading: ', trim(sfile)
      STOP
    ENDIF
    offset = szreal
    CALL MPI_FILE_READ_AT_ALL(fh, offset,rvar,1,GC_REAL,this%istatus_,this%ierr_) ! time
    offset = 2*szreal
    this%nparts_ = 0
    nb = 0
    nr = this%maxparts_/this%nprocs_
    DO WHILE ((this%ierr_.EQ.MPI_SUCCESS) .AND. (nb.LT.this%maxparts_))
      nr = MIN(nr, this%maxparts_-nb)
      CALL MPI_FILE_READ_AT_ALL(fh,offset,this%ptmp0_,3*nr,GC_REAL,this%istatus_,this%ierr_) ! PDB
      offset = offset + 3*nr*szreal
      DO j = 1,nr
        IF (this%wrtunit_ .EQ. 1) THEN
          this%ptmp0_(3,j) = this%ptmp0_(3,j)*this%invdel_(3)
        END IF
        IF ((this%ptmp0_(3,j).GE.this%lxbnds_(3,1)).AND.(this%ptmp0_(3,j).LT.this%lxbnds_(3,2))) THEN
          IF (this%nparts_.GE.this%partbuff_) THEN
            this%partbuff_ = this%partbuff_ + this%partchunksize_
            CALL this%ResizeArrays(this%partbuff_,.true.)
          END IF
          this%nparts_ = this%nparts_+1
          this%id_(this%nparts_) = j+nb-1
        END IF
      END DO
      nb = nb + nr
    END DO
    CALL MPI_FILE_CLOSE(fh,this%ierr_)
  END SUBROUTINE binary_read_id_co

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : binary_read_pdb_co
  !!  DESCRIPTION: Does read of binary Lagrangian particle data from file, 
  !!               collectively.
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance
  !!    iunit   : unit number
  !!    sfile   : fully resolved file name 
  !!    pdb     : part. d.b.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE binary_read_pdb_co(this,iunit,sfile,time,pdb)
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT)        :: this
    REAL(KIND=GP)                             :: rvar,time
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:,:):: pdb
    INTEGER,INTENT(IN)                        :: iunit
    INTEGER                                   :: i,j,nerr,szreal,nr,nb
    INTEGER(kind=MPI_OFFSET_KIND)             :: offset
    TYPE(MPI_File)                            :: fh
    CHARACTER(len=*),INTENT   (IN)            :: sfile

    CALL MPI_TYPE_SIZE(GC_REAL,szreal,this%ierr_)
    CALL MPI_FILE_OPEN(this%comm_,trim(sfile),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,this%ierr_)
    IF ( this%ierr_ .NE. MPI_SUCCESS ) THEN
      CALL MPI_ERROR_STRING(this%ierr_, this%serr_, nerr,ierr);
      WRITE(*,*) 'binary_read_pdb_co: Error reading opening : ', trim(sfile),& 
                trim(this%serr_)
      STOP
    ENDIF
  
    ! Must read part. data from correct spot in file:
    offset = 0
    CALL MPI_FILE_READ_AT_ALL(fh,offset,rvar,1,GC_REAL,this%istatus_,this%ierr_)    !  no.parts
    IF ( int(rvar).NE.this%maxparts_ ) THEN
      WRITE(*,*) 'binary_read_pdb_co: Attempt to read incorrect number of particles: required:', &
                  this%maxparts_,' no. read: ',int(rvar)
      WRITE(*,*) 'binary_read_pdb_co: Error reading: ', trim(sfile)
      STOP
    ENDIF
    offset = szreal
    CALL MPI_FILE_READ_AT_ALL(fh, offset,rvar,1,GC_REAL,this%istatus_,this%ierr_) ! time
    offset = 2*szreal
    IF (this%iexchtype_.EQ.GPEXCHTYPE_NN) THEN
      i = 1
      nb = 0
      nr = this%maxparts_/this%nprocs_
      DO WHILE ((this%ierr_.EQ.MPI_SUCCESS) .AND. (nb.LT.this%maxparts_))
        nr = MIN(nr, this%maxparts_-nb)
        CALL MPI_FILE_READ_AT_ALL(fh,offset,this%gptmp0_,3*nr,GC_REAL,this%istatus_,this%ierr_) ! PDB
        offset = offset + 3*nr*szreal
        DO j = 1,nr
          IF ((i.LE.this%nparts_).AND.(this%id_(i).EQ.(j+nb-1))) THEN
            pdb(:,i) = this%gptmp0_(:,j)
            i = i + 1
          END IF
        END DO
        nb = nb + nr
      END DO
    ELSE IF (this%iexchtype_.EQ.GPEXCHTYPE_VDB) THEN
      CALL MPI_FILE_READ_AT_ALL(fh,offset,pdb,3*this%maxparts_,GC_REAL,this%istatus_,this%ierr_) ! PDB
    END IF
    CALL MPI_FILE_CLOSE(fh,this%ierr_)
  END SUBROUTINE binary_read_pdb_co


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : binary_read_pdb_t0
  !!  DESCRIPTION: Does binary read of Lagrangian position d.b. from file
  !!               only from MPI task 0, and broadcast to all other tasks.
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance
  !!    iunit   : unit number
  !!    sfile   : fully resolved file name
  !!    time    : real time
  !!    pdb     : part. d.b. in array
  !!    stg     : stage of reading (if True, only determine ids from 
  !!                                file and resize if necessary)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE binary_read_pdb_t0(this,iunit,sfile,time,pdb,stg)
    USE fprecision
    USE commtypes
    USE mpivars
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT)         :: this
    REAL(KIND=GP),INTENT(INOUT)                :: time
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:,:) :: pdb
    LOGICAL      ,INTENT(IN), OPTIONAL         :: stg
    REAL(KIND=GP)                              :: fnt
    INTEGER      ,INTENT   (IN)                :: iunit
    INTEGER                                    :: j
    CHARACTER(len=*),INTENT(IN)                :: sfile
    LOGICAL                                    :: calc_ids

    ! Read global VDB, with time header, indexed only
    ! by time index: dir/spref.TTT.lag:
    IF ( this%myrank_.EQ.0 ) THEN
      OPEN(iunit,file=trim(sfile),status='old',access='stream', &
           form='unformatted',iostat=this%ierr_)
      IF ( this%ierr_.NE.0 ) THEN
        WRITE(*,*)'binary_read_pdb_t0: could not open file for reading: ',&
        trim(sfile)
        STOP
      ENDIF

      REWIND(iunit)
      READ(iunit) fnt
      READ(iunit) time
      IF ( int(fnt).NE.this%maxparts_ ) THEN
        WRITE(*,*)this%myrank_, &
          ': binary_read_pdb_t0: particle inconsistency: no. required=',&
          this%maxparts_,' no. found=',int(fnt), &
          ' file=',trim(sfile)
        STOP
      ENDIF
      READ(iunit) pdb
      CLOSE(iunit)
    ENDIF

    IF (this%iexchtype_.EQ.GPEXCHTYPE_VDB) THEN
      CALL MPI_BCAST(pdb,3*this%maxparts_,GC_REAL,0,this%comm_,this%ierr_)
      IF ( this%ierr_.NE.MPI_SUCCESS ) THEN
        WRITE(*,*)this%myrank_, ': binary_read_pdb_t0: Broadcast failed: file=',&
        trim(sfile)
    ENDIF
    ELSE IF (this%iexchtype_.EQ.GPEXCHTYPE_NN) THEN
      IF (this%myrank_.EQ.0) THEN
!$omp parallel do
        DO j = 1,this%maxparts_
          this%id_(j) = j-1
        END DO
      END IF
      IF (PRESENT(stg)) THEN
        calc_ids = stg
      ELSE
        calc_ids = .false.
      END IF
      IF (calc_ids) THEN
        IF (this%myrank_.EQ.0) this%nparts_ = this%maxparts_
        IF (this%wrtunit_ .EQ.1) pdb(3,:) = pdb(3,:)*this%invdel_(3)
        CALL this%gpcomm_%IdentifyTaskV(this%id_,pdb(3,:),this%nparts_,this%tmpint_)
        IF (this%nparts_.GT.this%partbuff_) THEN
          PRINT *, 'Rank', this%myrank_, 'resizing: nparts=', this%nparts_,  &
                   ' | partbuff=', this%partbuff_, ' --> ', &
                   (1+this%nparts_/this%partchunksize_)*this%partchunksize_
          this%partbuff_ = (1+this%nparts_/this%partchunksize_)*this%partchunksize_
          CALL this%ResizeArrays(this%partbuff_,.true.)
        END IF
      ELSE
        CALL this%gpcomm_%PartScatterV(this%id_,pdb(1,:),pdb(2,:),pdb(3,:),this%nparts_,this%tmpint_)
      END IF
    END IF
  END SUBROUTINE binary_read_pdb_t0


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : ascii_read_pdb
  !!  DESCRIPTION: Does ASCII read of Lagrangian position d.b. from file.
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance
  !!    iunit   : unit number
  !!    sfile   : fully resolved file name
  !!    time    : real time
  !!    pdb     : part. d.b. in (4,maxparts) array
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ascii_read_pdb(this,iunit,sfile,time,pdb)
    USE fprecision
    USE commtypes
    USE mpivars
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT) :: this
    REAL(KIND=GP),INTENT(INOUT)        :: time
    REAL(KIND=GP),INTENT(INOUT)        :: pdb(3,this%maxparts_)
    INTEGER      ,INTENT   (IN)        :: iunit
    INTEGER                            :: j,nt
    CHARACTER(len=*),INTENT(IN)        :: sfile

    ! Read global VDB, with time header, indexed only
    ! by time index: dir/spref.TTT.txt:
    IF ( this%myrank_.EQ.0 ) THEN
      OPEN(iunit,file=trim(sfile),status='old',form='formatted',iostat=this%ierr_)
      IF ( this%ierr_.NE.0 ) THEN
        WRITE(*,*)'ascii_read_pdb: could not open file for reading: ', &
        trim(sfile)
        STOP
      ENDIF
      READ(iunit,*,iostat=this%ierr_) nt
      READ(iunit,*,iostat=this%ierr_) time
      IF ( nt.LT.this%maxparts_ ) THEN
        WRITE(*,*)this%myrank_, &
          ': ascii_read_pdb: particle inconsistency: no. required=',   &
          this%maxparts_,' no. found=',nt, &
          ' file=',trim(sfile)
        STOP
      ENDIF
      DO j = 1, this%maxparts_
        READ(iunit,*,iostat=this%ierr_) pdb(1,j),pdb(2,j),pdb(3,j)
  600   FORMAT(3(E23.15,1X))
      ENDDO
      CLOSE(iunit)
    ENDIF
    CALL MPI_BCAST(pdb,3*this%maxparts_,GC_REAL,0,this%comm_,this%ierr_)
    IF ( this%ierr_.NE.MPI_SUCCESS ) THEN
        WRITE(*,*)this%myrank_, ': ascii_read_pdb: Broadcast failed: file=', &
        trim(sfile)
    ENDIF
  END SUBROUTINE ascii_read_pdb


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : io_readvec
  !!  DESCRIPTION: Reads a 3-component Lagrangian vector field from
  !!               file (same binary/ASCII format as xlg) and scatters
  !!               the three components into pstate(start_index  ) (+1 +2).
  !!               Must be called AFTER io_read has populated id_,
  !!               nparts_, and (for VDB) vdb_, so that the correct
  !!               local-to-global particle mapping is already known.
  !!               This is the main entry point for both binary and
  !!               ASCII reads.
  !!  ARGUMENTS  :
  !!    this        : 'this' class instance
  !!    pstate      : particle state array (read positions must already
  !!                  be loaded; provides context for NN id_ mapping)
  !!    iunit       : Fortran unit number to use for sequential I/O
  !!    dir         : input directory
  !!    spref       : filename prefix (e.g. 'vip')
  !!    nmb         : time index string; if len_trim(nmb)==0, spref is
  !!                  treated as a fully-resolved filename and dir is
  !!                  ignored (same convention as io_read)
  !!    start_index : first pstate index to fill (e.g. this%VELOCITY)
  !!    opiotype    : optional; overrides iouttype_ if present
  !!    opbcoll     : optional; overrides bcollective_ if present
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE io_readvec(this, pstate, iunit, dir, spref, nmb, &
                        start_index, opiotype, opbcoll)
    USE gpstate_mod
    USE fprecision
    USE commtypes
    USE mpivars
!$  USE threads
    IMPLICIT NONE
    CLASS(ParticleBase)    ,INTENT(INOUT)          :: this
    TYPE  (GPStateComp)    ,INTENT(INOUT)          :: pstate(:)
    INTEGER,INTENT(IN)                             :: iunit
    INTEGER,INTENT(IN)                             :: start_index
    INTEGER,INTENT(IN),OPTIONAL                    :: opiotype, opbcoll
    CHARACTER(len=128),INTENT(IN)                  :: dir
    CHARACTER(len=*)  ,INTENT(IN)                  :: nmb
    CHARACTER(len=*)  ,INTENT(IN)                  :: spref
    REAL(KIND=GP)                                  :: time
    INTEGER                                        :: bcoll, iotype, j, ng
    INTEGER                                        :: ht

    ! Resolve I/O mode overrides
    IF ( present(opiotype) ) THEN
      iotype = opiotype
    ELSE
      iotype = this%iouttype_
    ENDIF
    IF ( present(opbcoll) ) THEN
      bcoll = opbcoll
    ELSE
      bcoll = this%bcollective_
    ENDIF

    ! Read global vector data into gptmp0_(3, maxparts_).
    ! We reuse the same low-level readers as io_read; they expect a
    ! (3, maxparts_) scratch array and populate it in global particle-
    ! id order. For the NN collective path we read into gptmp0_ via
    ! binary_read_pdb_co (which already handles the NN id-based
    ! scatter into a local array); for all other paths we call the
    ! same t0/ASCII readers used by io_read.
    CALL GTStart(this%htimers_(GPTIME_GPREAD))
    IF ( iotype .EQ. 0 ) THEN   ! Binary files
      IF ( bcoll .EQ. 1 ) THEN  ! Collective binary
        IF ( this%iexchtype_ .EQ. GPEXCHTYPE_NN ) THEN
          ! binary_read_pdb_co handles the NN case internally: it
          ! reads chunks collectively and keeps only the records
          ! whose global index matches a local id_, storing them
          ! compactly in gptmp0_(:, 1:nparts_).
          IF (len_trim(nmb).gt.0) THEN
            CALL this%binary_read_pdb_co(iunit,                           &
                 trim(dir) // '/' // trim(spref) // '.' // nmb // '.lag', &
                 time, this%gptmp0_)
          ELSE
            CALL this%binary_read_pdb_co(iunit, trim(spref),              &
                 time, this%gptmp0_)
          ENDIF
          ! gptmp0_(:, j) now holds the vector for local particle j.
          ! We copy directly into pstate slots.
!$omp parallel do
          DO j = 1, this%nparts_
            pstate(start_index  )%rcomp(j) = this%gptmp0_(1,j)
            pstate(start_index+1)%rcomp(j) = this%gptmp0_(2,j)
            pstate(start_index+2)%rcomp(j) = this%gptmp0_(3,j)
          ENDDO
        ELSE  ! GPEXCHTYPE_VDB, collective binary
          IF (len_trim(nmb).gt.0) THEN
            CALL this%binary_read_pdb_co(iunit,                           &
                 trim(dir) // '/' // trim(spref) // '.' // nmb // '.lag', &
                 time, this%gptmp0_)
          ELSE
            CALL this%binary_read_pdb_co(iunit, trim(spref),              &
                 time, this%gptmp0_)
          ENDIF
          ! gptmp0_ holds the full global vector; scatter via vdb_.
          CALL this%CopyLocalWrk(pstate(start_index  )%rcomp,             &
                                 pstate(start_index+1)%rcomp,             &
                                 pstate(start_index+2)%rcomp,             &
                                 this%vdb_, this%gptmp0_, this%maxparts_)
        ENDIF  ! iexchtype
      ELSE  ! Non-collective (task-0) binary
        IF (len_trim(nmb).gt.0) THEN
          CALL this%binary_read_pdb_t0(iunit,                             &
               trim(dir) // '/' // trim(spref) // '.' // nmb // '.lag',   &
               time, this%gptmp0_)
        ELSE
          CALL this%binary_read_pdb_t0(iunit, trim(spref),                &
               time, this%gptmp0_)
        ENDIF
        ! After binary_read_pdb_t0: for VDB the global data is
        ! broadcast and sits in gptmp0_; for NN the data has been
        ! scatter-communicated, so gptmp0_(:, 1:nparts_) is local.
        CALL this%io_readvec_scatter_(pstate, start_index)
      ENDIF  ! bcoll
    ELSE  ! ASCII files
      IF (len_trim(nmb).gt.0) THEN
        CALL this%ascii_read_pdb(iunit,                                   &
             trim(dir) // '/' // trim(spref) // '.' // nmb // '.txt',     &
             time, this%gptmp0_)
      ELSE
        CALL this%ascii_read_pdb(iunit, trim(spref),                      &
             time, this%gptmp0_)
      ENDIF
      ! ascii_read_pdb always broadcasts to all ranks, so gptmp0_
      ! holds the full global array on every rank.
      CALL this%io_readvec_scatter_(pstate, start_index)
    ENDIF  ! iotype
    CALL GTAcc(this%htimers_(GPTIME_GPREAD))
    ! Sanity check: global particle count must still equal maxparts_
    CALL MPI_ALLREDUCE(this%nparts_, ng, 1, MPI_INTEGER,                  &
                       MPI_SUM, this%comm_, this%ierr_)
    IF ( this%myrank_.EQ.0 .AND. ng.NE.this%maxparts_ ) THEN
      WRITE(*,*) 'io_readvec: inconsistent d.b. after read of ', trim(spref), &
                 ': expected: ', this%maxparts_, '; found: ', ng
      STOP
    ENDIF
  END SUBROUTINE io_readvec


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : io_readvec_scatter_ (private helper)
  !!  DESCRIPTION: After a non-collective (t0 or ASCII) read has
  !!               filled gptmp0_(3, maxparts_) with the global
  !!               vector data (broadcast to all ranks for VDB,
  !!               scatter-communicated for NN), copy the locally-
  !!               owned records into the three pstate slots.
  !!               For VDB: gptmp0_ is a full global array; we walk
  !!               vdb_ to identify which particles this rank owns,
  !!               exactly as CopyLocalWrk does, but into pstate.
  !!               For NN:  binary_read_pdb_t0 already called
  !!               PartScatterV internally, so gptmp0_(:,1:nparts_)
  !!               already contains only local data in compact form.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE io_readvec_scatter_(this, pstate, start_index)
    USE gpstate_mod
    USE fprecision
    USE commtypes
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT) :: this
    TYPE  (GPStateComp) ,INTENT(INOUT) :: pstate(:)
    INTEGER             ,INTENT(IN)    :: start_index
    INTEGER                            :: j

    IF ( this%iexchtype_ .EQ. GPEXCHTYPE_VDB ) THEN
      ! gptmp0_ is a full (3, maxparts_) global array on every rank.
      ! Use vdb_ (which records the position of every particle in
      ! global id order) to decide which belong to this z-slab.
      CALL this%CopyLocalWrk(pstate(start_index  )%rcomp,              &
                             pstate(start_index+1)%rcomp,              &
                             pstate(start_index+2)%rcomp,              &
                             this%vdb_, this%gptmp0_, this%maxparts_)
    ELSE  ! GPEXCHTYPE_NN
      ! binary_read_pdb_t0 called PartScatterV for the position data, but
      ! that scatter applied to ptmp0_, not gptmp0_.  However, for the
      ! NN/non-collective path, task 0 holds the entire file and PartScatterV
      ! distributes each particle to the rank that owns its z-slab.  For a 
      ! *vector* file we followed the same binary_read_pdb_t0 call, which 
      ! internally called PartScatterV on our gptmp0_ array. So
      ! gptmp0_(:,1:nparts_) now contains only the locally-owned records.
!$omp parallel do
      DO j = 1, this%nparts_
        pstate(start_index  )%rcomp(j) = this%gptmp0_(1,j)
        pstate(start_index+1)%rcomp(j) = this%gptmp0_(2,j)
        pstate(start_index+2)%rcomp(j) = this%gptmp0_(3,j)
      ENDDO

    ENDIF
  END SUBROUTINE io_readvec_scatter_  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : EulerToLag
  !!  DESCRIPTION: Computes the Eulerian to Lagrangian
  !!               transformation by interpolating Eulerian field
  !!               evar to current position of Lagrangian paricles in 
  !!               d.b. Array lag must be large enough to accommodate 
  !!               max. no. particles; no checking is done. Note
  !!               that 'evar' array must have local dimensions 
  !!               for a real array in GHOST (nx X ny X (kend-ksta+1)).
  !!               Global computation of spline or other interpolation 
  !!               operator will be done here.
  !!
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance (IN)
  !!    lag     : real array with the containing Lagrangian value
  !!              of evar field on output (OUT)
  !!    nl      : no. Lag. points in lag
  !!    evar    : Eulerian variable (IN)
  !!    doupdate: if true, do interp point update in interpolator; else don't
  !!    tmp1/2  : temp. arrays for interpolation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE EulerToLag(this,lag,nl,evar,doupdate,tmp1,tmp2)
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT)                     :: this
    INTEGER      ,INTENT   (IN)                            :: nl
    LOGICAL      ,INTENT   (IN)                            :: doupdate
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: evar,tmp1,tmp2
    REAL(KIND=GP),INTENT(INOUT),DIMENSION             (nl) :: lag
    INTEGER                                                :: j

    IF ( doupdate ) THEN
      CALL GTStart(this%htimers_(GPTIME_PUPDATE))
      CALL this%intop_%PartUpdate3D(this%px_,this%py_,this%pz_,this%nparts_)
      CALL GTAcc(this%htimers_(GPTIME_PUPDATE))
    ENDIF
    CALL GTStart(this%htimers_(GPTIME_SPLINE))
    CALL this%intop_%CompSpline3D(evar,tmp1,tmp2)
    CALL GTAcc(this%htimers_(GPTIME_SPLINE))

    CALL GTStart(this%htimers_(GPTIME_INTERP))
    CALL this%intop_%DoInterp3D(lag,nl)
    CALL GTAcc(this%htimers_(GPTIME_INTERP))
  END SUBROUTINE EulerToLag

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : MakePeriodicP
  !!  DESCRIPTION: Enforces periodic b.c.'s on particles in pdb
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance (IN)
  !!    px      : particle x loc. d.b.
  !!    py      : particle y loc. d.b.
  !!    pz      : particle z loc. d.b.
  !!    npdb    : no. particles in pdb
  !!    idir    : first three bits provide which directions to periodize.
  !!              So, 1==>x, 2==>y, 4==>z; 3==>x & y; 7==>x & y & z
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE MakePeriodicP(this,px,py,pz,npdb,idir)
    USE fprecision
    USE commtypes
    USE mpivars
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT)               :: this
    INTEGER,INTENT(IN)                               :: idir,npdb
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)         :: px,py,pz
    INTEGER                                          :: j

! The old form of the wrap, modulo(p+2*L,L), also perturbs by one ULP
! the low bits of most particles that do NOT wrap (adding 2L and
! reducing is not bit-exact), and it did so at every stage of every
! step. Besides the noise this injected in the trajectories, the VDB
! exchange periodizes all three directions while the NN exchange only
! periodizes z conditionally, so the two exchange modes slowly drifted
! apart from these seeds. The conditional form below, the same one
! MakePeriodicZ uses, leaves particles inside the box bit-untouched;
! it assumes the positions are within one box length of the domain,
! which holds at every call site since the callers periodize right
! after a single stage displacement (the NN exchange already requires
! motion of less than one zone per step).
! One parallel region for the three directions; every thread sees the
! same idir, so all of them take the same branches of the IFs.
!$omp parallel
    IF ( btest(idir,0) ) THEN
!$omp do
      DO j = 1, npdb
        IF ( px(j).LT.0 ) THEN
          px(j) = px(j) + this%gext_(1)
          ! p+L can round up to exactly L when p is a tiny negative
          ! number, and L is outside [0,L): fold it to 0. The subtract
          ! branch below is exact (Sterbenz) and needs no guard.
          IF ( px(j).GE.this%gext_(1) ) px(j) = 0.0_GP
        ELSE IF ( px(j).GE.this%gext_(1) ) THEN
          px(j) = px(j) - this%gext_(1)
        ENDIF
      ENDDO
    ENDIF
    IF ( btest(idir,1) ) THEN
!$omp do
      DO j = 1, npdb
        IF ( py(j).LT.0 ) THEN
          py(j) = py(j) + this%gext_(2)
          ! p+L can round up to exactly L when p is a tiny negative
          ! number, and L is outside [0,L): fold it to 0. The subtract
          ! branch below is exact (Sterbenz) and needs no guard.
          IF ( py(j).GE.this%gext_(2) ) py(j) = 0.0_GP
        ELSE IF ( py(j).GE.this%gext_(2) ) THEN
          py(j) = py(j) - this%gext_(2)
        ENDIF
      ENDDO
    ENDIF
    IF ( btest(idir,2) ) THEN
!$omp do
      DO j = 1, npdb
        IF ( pz(j).LT.0 ) THEN
          pz(j) = pz(j) + this%gext_(3)
          ! p+L can round up to exactly L when p is a tiny negative
          ! number, and L is outside [0,L): fold it to 0. The subtract
          ! branch below is exact (Sterbenz) and needs no guard.
          IF ( pz(j).GE.this%gext_(3) ) pz(j) = 0.0_GP
        ELSE IF ( pz(j).GE.this%gext_(3) ) THEN
          pz(j) = pz(j) - this%gext_(3)
        ENDIF
      ENDDO
    ENDIF
!$omp end parallel
  END SUBROUTINE MakePeriodicP

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : MakePeriodicZ
  !!  DESCRIPTION: Enforces periodic b.c.'s on particles in pdb
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance (IN)
  !!    pz      : particle z loc. d.b.
  !!    tpz     : particle z stored loc. d.b.
  !!    npdb    : no. particles in pdb
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE MakePeriodicZ(this,pz,tpz,npdb)
    USE fprecision
    USE commtypes
    USE mpivars
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT)               :: this
    INTEGER,INTENT(IN)                               :: npdb
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(npdb)      :: pz,tpz
    INTEGER                                          :: j
!$omp parallel do 
    ! Unlike MakePeriodicP, a result of p+L that rounds up to exactly L
    ! must NOT be folded to zero here. This routine wraps particles the
    ! NN exchange has just delivered to the opposite end of the domain,
    ! and the particle stays on the receiving task: on the task that
    ! owns the top of the box the seam must be represented as z=L,
    ! which lies inside its ghost-extended interpolation range, while
    ! z=0 does not. The VDB exchange, which re-bins ownership from the
    ! coordinate value, needs the opposite convention and gets it from
    ! MakePeriodicP.
    DO j = 1,npdb
       IF (pz(j).LT.0) THEN
          pz(j)  =  pz(j) + this%gext_(3)
          tpz(j) = tpz(j) + this%gext_(3)
       ELSE IF (pz(j).GE.this%gext_(3)) THEN
          pz(j)  =  pz(j) - this%gext_(3)
          tpz(j) = tpz(j) - this%gext_(3)
       ENDIF
    ENDDO
  END SUBROUTINE MakePeriodicZ

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : MakePeriodicExt
  !!  DESCRIPTION: Enforces periodic b.c.'s on extended field
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance (IN)
  !!    v       : real field on extended grid 
  !!    nx,ny   : global size of field in x-y (including ghost zones)
  !!    kb,ke   : starting, ending z-indices of slab (including ghost zones)
  !!    nc      : index in x and y (and z) s.t. f(nc+1,:,:) = f(nx-nc,:,:), etc.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE MakePeriodicExt(this,v,nx,ny,kb,ke,nc)
    USE fprecision
    USE commtypes
!$  USE threads
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT)               :: this
    INTEGER      ,INTENT(IN)                         :: nc,nx,ny,kb,ke
    REAL(KIND=GP),INTENT(INOUT)                      :: v(nx,ny,kb:ke)
    INTEGER                                          :: i,j,k
    ! Recall: nx, ny are the dimensions _including_ ghost zones:
    !
    ! Periodicity s.t.:
    !   | | [ | | | | ] | |
    !   a b       a b 
!$omp parallel do if (ke-kb.ge.nth) private (i,j,k)
    DO k = kb,ke 
!$omp parallel do if (ke-kb.lt.nth) private (i,j)
      DO j = 1,ny
        DO i = 1,nc
          v(i,j,k) = v(nx-nc+i,j,k)
        ENDDO
        DO i = nx-nc+1,nx
          v(i,j,k) = v(2*nc+i-nx,j,k)
        ENDDO
      ENDDO
      DO i = 1,nx
        DO j = 1,nc
          v(i,j,k) = v(i,nx-nc+j,k)
        ENDDO
        DO j = ny-nc+1,ny
          v(i,j,k) = v(i,2*nc+j-nx,k)
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE MakePeriodicExt

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : Part_Delete
  !!  DESCRIPTION: Removes from PDB NULL particles, concatenates list,
  !!               and sets new number of particles
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance (IN)
  !!    id      : part ids
  !!    px,py pz: part. d.b.
  !!    npdb    : no. parts. in pdb
  !!    nnew    : no. non-NULL particles (set in GPartComm class)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Part_Delete(this,id,px,py,pz,npdb,nnew)
    USE fprecision
    USE commtypes
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT)         :: this
    INTEGER      ,INTENT   (IN)                :: npdb
    INTEGER      ,INTENT(INOUT),DIMENSION(npdb):: id
    INTEGER      ,INTENT  (OUT)                :: nnew
    INTEGER                                    :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(npdb):: px,py,pz
    nnew = 0
!$omp parallel do 
    DO i = 1, npdb
       IF ( this%id_(i) .NE. GPNULL ) nnew = nnew + 1
    ENDDO
    j    = 1
    DO i = 1, nnew
      DO WHILE ( j.LE.npdb .AND. id(j).EQ.GPNULL )
        j = j + 1
      ENDDO
      IF ( j.LE.npdb .AND. j.NE.i ) THEN
        id(i) = id(j); id(j) = GPNULL
        px(i) = px(j)
        py(i) = py(j)
        pz(i) = pz(j)
      ENDIF
    ENDDO
  END SUBROUTINE Part_Delete

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : SelectLocalIds
  !!  DESCRIPTION: Builds the list of the particles whose z position
  !!               in gvdb lies in [zlo,zhi), in ascending id order.
  !!               The ascending order is what every consumer of the
  !!               local arrays assumes (CopyLocalWrk and the io
  !!               scatter helpers walk the VDB in ascending order),
  !!               so it must be deterministic: the previous version
  !!               of this selection, threaded with the counter in a
  !!               critical section, filled the arrays in whatever
  !!               order the threads arrived, which mislabeled
  !!               particles. The selection here is done in two
  !!               passes, so each thread counts the matches of its
  !!               contiguous chunk and then writes them at its
  !!               offset: the result is identical to a sequential
  !!               ascending scan for any number of threads (and for
  !!               builds without OpenMP), and it is also faster than
  !!               both the sequential scan and the old critical
  !!               section version.
  !!  ARGUMENTS  :
  !!    zlo,zhi : z-slab bounds; a particle is local if
  !!              zlo <= gvdb(3,j) < zhi
  !!    gvdb    : global VDB with particle position records
  !!    ngvdb   : no. records in gvdb
  !!    id      : on output, id(1:nl) holds the 0-based ids of the
  !!              local particles, ascending
  !!    nl      : no. local particles found
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SelectLocalIds(zlo,zhi,gvdb,ngvdb,id,nl)
    USE fprecision
!$  USE omp_lib
    IMPLICIT NONE
    REAL(KIND=GP),INTENT   (IN)                    :: zlo,zhi
    INTEGER      ,INTENT   (IN)                    :: ngvdb
    REAL(KIND=GP),INTENT   (IN),DIMENSION(3,ngvdb) :: gvdb
    INTEGER      ,INTENT(INOUT),DIMENSION(*)       :: id
    INTEGER      ,INTENT  (OUT)                    :: nl
    INTEGER                                        :: j,k,t,lo,hi,chunk,ntmax
    INTEGER                                        :: nthr
    INTEGER      ,ALLOCATABLE  ,DIMENSION(:)       :: cnt,base

    ntmax = 1
!$  ntmax = omp_get_max_threads()
    ALLOCATE(cnt(0:ntmax-1),base(0:ntmax-1))
    nthr = 1
!$omp parallel private(t,j,k,lo,hi,chunk) shared(cnt,base,nthr)
!$omp single
!$  nthr = omp_get_num_threads()
!$omp end single
    t = 0
!$  t = omp_get_thread_num()
    chunk = (ngvdb+nthr-1)/nthr
    lo = t*chunk + 1
    hi = MIN(ngvdb,(t+1)*chunk)
    ! Pass 1: each thread counts the matches in its chunk
    k = 0
    DO j = lo, hi
      IF ( gvdb(3,j).GE.zlo .AND. gvdb(3,j).LT.zhi ) k = k + 1
    ENDDO
    cnt(t) = k
!$omp barrier
!$omp single
    base(0) = 0
    DO j = 1, nthr-1
      base(j) = base(j-1) + cnt(j-1)
    ENDDO
!$omp end single
    ! Pass 2: each thread writes its matches at its offset
    k = base(t)
    DO j = lo, hi
      IF ( gvdb(3,j).GE.zlo .AND. gvdb(3,j).LT.zhi ) THEN
        k = k + 1
        id(k) = j - 1
      ENDIF
    ENDDO
!$omp end parallel
    nl = base(nthr-1) + cnt(nthr-1)
    DEALLOCATE(cnt,base)
  END SUBROUTINE SelectLocalIds


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : GetLocalWrk
  !!  DESCRIPTION: Removes from PDB NULL particles, concatenates list,
  !!               and sets new number of particles
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance (IN)
  !!    id      : part ids, returned if gfill not specified
  !!    lx,ly,lz: local part. d.b. vectors
  !!    nl      : no. parts. in local pdb. If gfill specified, this is read
  !!              as the local no. particles.
  !!    gvdb    : global VDB containing part. position records. Location
  !!              gives particle id.
  !!    ngvdb   : no. records in global VDB
  !!    gfill   : if specified, the gvdb will be used to locate the 
  !!              particles this task owns, and will return the correct
  !!              local arrays, lx, ly, lz, from the global d.b. gfill,
  !!              using id array as indirection.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE GetLocalWrk(this,id,lx,ly,lz,nl,gvdb,ngvdb,gfill)
    USE fprecision
    USE commtypes
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT)                     :: this
    INTEGER      ,INTENT(INOUT)                            :: nl
    INTEGER      ,INTENT(INOUT),DIMENSION(this%maxparts_)  :: id
    INTEGER      ,INTENT   (IN)                            :: ngvdb
    INTEGER                                                :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(this%maxparts_)  :: lx,ly,lz
    REAL(KIND=GP),INTENT   (IN),DIMENSION(3,ngvdb)         :: gvdb
    REAL(KIND=GP),INTENT   (IN),DIMENSION(3,ngvdb),OPTIONAL:: gfill
    IF ( .NOT.present(gfill) ) THEN
      id = GPNULL
      ! Deterministic (ascending id) selection of the local particles
      CALL SelectLocalIds(this%lxbnds_(3,1),this%lxbnds_(3,2),gvdb,ngvdb,id,nl)
!$omp parallel do private (i)
      DO j = 1, nl
        i = id(j) + 1
        lx (j) = gvdb(1,i)
        ly (j) = gvdb(2,i)
        lz (j) = gvdb(3,i)
      ENDDO
    ELSE
!$omp parallel do
      DO j = 1, nl
        lx (j) = gfill(1,id(j)+1)
        ly (j) = gfill(2,id(j)+1)
        lz (j) = gfill(3,id(j)+1)
      ENDDO
    ENDIF
  END SUBROUTINE GetLocalWrk

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : GetLocalWrk_aux
  !!  DESCRIPTION: Removes from PDB NULL particles, concatenates list,
  !!               and sets new number of particles. This auxiliary 
  !!               subroutines also updates arrays used during the 
  !!               intermediate steps of the RK solver, and is needed 
  !!               if local work is recomputed in the midde of a RK 
  !!               iteration.
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance (IN)
  !!    id      : part ids
  !!    lx,ly,lz: local part. d.b. vectors
  !!    tx,ty,tz: local initial part. d.b. vectors
  !!    nl      : no. parts. in local pdb
  !!    gvdb    : global VDB containing part. position records. Location
  !!              gives particle id.
  !!    gtmp    : global VDB containing part. position records at the
  !!              beginning of the RK loop. Location gives particle id.
  !!    ngvdb   : no. records in global VDB
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE GetLocalWrk_aux(this,id,lx,ly,lz,tx,ty,tz,nl,gvdb,gtmp,ngvdb)
    USE fprecision
    USE commtypes
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT)                    :: this
    INTEGER      ,INTENT(INOUT)                           :: nl
    INTEGER      ,INTENT(INOUT),DIMENSION(this%maxparts_) :: id
    INTEGER      ,INTENT   (IN)                           :: ngvdb
    INTEGER                                               :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(this%maxparts_) :: lx,ly,lz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(this%maxparts_) :: tx,ty,tz
    REAL(KIND=GP),INTENT   (IN),DIMENSION(3,ngvdb)        :: gvdb,gtmp
    id = GPNULL
    ! Deterministic (ascending id) selection of the local particles
    CALL SelectLocalIds(this%lxbnds_(3,1),this%lxbnds_(3,2),gvdb,ngvdb,id,nl)
!$omp parallel do private (i)
    DO j = 1, nl
      i = id(j) + 1
      lx (j) = gvdb(1,i)
      ly (j) = gvdb(2,i)
      lz (j) = gvdb(3,i)
      tx (j) = gtmp(1,i)
      ty (j) = gtmp(2,i)
      tz (j) = gtmp(3,i)
    ENDDO
  END SUBROUTINE GetLocalWrk_aux

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : CopyLocalWrk
  !!  DESCRIPTION: Updates records of the VDB.
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance (IN)
  !!    lx,ly,lz: local part. d.b. vectors
  !!    gvdb    : global VDB containing part. position records. Location
  !!              gives particle id.
  !!    vgvdb   : global VDB containing part. property records
  !!              (can be velocity or anything associated to the particle).
  !!    ngvdb   : no. records in global VDB
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE CopyLocalWrk(this,lx,ly,lz,gvdb,vgvdb,ngvdb)
    USE fprecision
    USE commtypes
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT)                    :: this
    INTEGER      ,INTENT   (IN)                           :: ngvdb
    INTEGER                                               :: i,j,nll
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(this%maxparts_) :: lx,ly,lz
    REAL(KIND=GP),INTENT   (IN),DIMENSION(3,ngvdb)        :: gvdb
    REAL(KIND=GP),INTENT   (IN),DIMENSION(3,ngvdb)        :: vgvdb
    nll = 0
    DO j = 1, ngvdb
      IF ( gvdb(3,j).GE.this%lxbnds_(3,1) .AND. gvdb(3,j).LT.this%lxbnds_(3,2) ) THEN 
        nll = nll + 1
        lx (nll) = vgvdb(1,j)
        ly (nll) = vgvdb(2,j)
        lz (nll) = vgvdb(3,j)
      ENDIF
    ENDDO
  END SUBROUTINE CopyLocalWrk

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : GetVDB
  !!  DESCRIPTION: Gets particle d.b.
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance (IN)
  !!    pdb     : part pdb, of size (3,npdb)
  !!    npdb    : size of pdb array (2nd dimension); must be >= maxparts_
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE GetVDB(this,pdb,npdb)
    USE fprecision
    USE commtypes

    IMPLICIT NONE 
    CLASS(ParticleBase) ,INTENT(INOUT)            :: this 
    INTEGER      ,INTENT   (IN)                   :: npdb
    INTEGER                                       :: j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(3,npdb) :: pdb
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_NN ) THEN
      IF ( .NOT.this%PartNumConsistent(this%nparts_) ) THEN
          IF ( this%myrank_.eq.0 ) THEN
            WRITE(*,*) 'GetVDB: Inconsistent particle count'
            STOP
        ENDIF
      ENDIF
      CALL this%gpcomm_%VDBSynch(pdb,this%maxparts_,this%id_, &
           this%px_,this%py_,this%pz_,this%nparts_,this%ptmp0_)
    ELSE
!$omp parallel do 
      DO j = 1, npdb
        pdb(1:3,j) = this%vdb_(1:3,j)
     ENDDO
    ENDIF
   END SUBROUTINE GetVDB

   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : GetVel
  !!  DESCRIPTION: Gets current particle velocities by doing 'synch'
  !!               of local velocities
  !!         
  !!  ARGUMENTS  :
  !!    this     : 'this' class instance (IN)
  !!    lvel     : part velocity array, of size (3,nparts)
  !!    nparts   : size of lvel array, must be >= maxparts_
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE GetVel(this,lvel,nparts)
    USE fprecision
    USE commtypes
    IMPLICIT NONE 
    CLASS(ParticleBase) ,INTENT(INOUT)             :: this 
    INTEGER      ,INTENT   (IN)                    :: nparts
    INTEGER                                        :: j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(3,nparts):: lvel
    IF ( .NOT.this%PartNumConsistent(this%nparts_) ) THEN
      IF ( this%myrank_.eq.0 ) THEN
        WRITE(*,*) 'GetVel: Inconsistent particle count'
        STOP
      ENDIF
    ENDIF
    CALL this%gpcomm_%VDBSynch(lvel,this%maxparts_,this%id_, &
         this%lvx_,this%lvy_,this%lvz_,this%nparts_,this%ptmp0_)
   END SUBROUTINE GetVel

   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : GetTime
  !!  DESCRIPTION: gets elapsed time from timer index itime
  !!         
  !!  ARGUMENTS  :
  !!    this     : 'this' class instance (IN)
  !!    itime    : 'GPTIME' parameter (above)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION FUNCTION GetTime(this,itime)
    USE fprecision
    USE commtypes
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT)             :: this
    INTEGER      ,INTENT   (IN)                    :: itime
    INTEGER                                        :: j
    IF ( itime.LT.GPTIME_STEP .OR. itime.GT.GPMAXTIMERS ) THEN
      WRITE(*,*)'GetTime: invalid time specification'
      STOP
    ENDIF
    GetTime = GTGetTime(this%htimers_(itime))
   END FUNCTION GetTime

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : GetNParts
  !!  DESCRIPTION: Gets no. particles on grid
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance (IN)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER FUNCTION GetNParts(this)
     CLASS(ParticleBase) ,INTENT(INOUT)            :: this 
     INTEGER                                       :: ngp    
     CALL MPI_ALLREDUCE(this%nparts_,ngp,1,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
     GetNParts = ngp
  END FUNCTION GetNParts

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : GetLoadBal
  !!  DESCRIPTION: Gets current load (im)balance measure
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance (IN)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL FUNCTION GetLoadBal(this)
     CLASS(ParticleBase) ,INTENT(INOUT)            :: this 
     REAL                                          :: rbal      
     INTEGER                                       :: gnmax,gnmin    
     CALL MPI_ALLREDUCE(this%nparts_,gnmin,1,MPI_INTEGER,MPI_MIN,this%comm_,this%ierr_)
     CALL MPI_ALLREDUCE(this%nparts_,gnmax,1,MPI_INTEGER,MPI_MAX,this%comm_,this%ierr_)
     rbal = real(gnmax) / (real(gnmin)+tiny(1.0))
     GetLoadBal = rbal
  END FUNCTION GetLoadBal

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : PartNumConsistent
  !!  DESCRIPTION: Checks that sum of local particle counts equals
  !!               maxparts_
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance (IN)
  !!    nlocal  : local part. count
  !!    gsum    : (optional) global sum of all nlocal
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  LOGICAL FUNCTION PartNumConsistent(this,nlocal,gsum)
     CLASS(ParticleBase) ,INTENT(INOUT)            :: this 
     REAL                                          :: rbal      
     INTEGER,INTENT(IN)                            :: nlocal
     INTEGER,INTENT(OUT),OPTIONAL                  :: gsum
     INTEGER                                       :: ng
     CALL MPI_ALLREDUCE(nlocal,ng,1,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
     IF ( present(gsum) ) gsum = ng
     PartNumConsistent = ng .EQ. this%maxparts_
  END FUNCTION PartNumConsistent

 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : R3toR3
  !!  DESCRIPTION: Copies input 3D real array to output 3D real array.
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance
  !!    vout    : result, returned; size standard in GHOST: (nx,ny,ksta:kend)
  !!    vin     : input array, size standard in GHOST
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE R3toR3(this, vout, vin)
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars
!$  USE threads
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT)                   :: this
    REAL(KIND=GP),INTENT(OUT),DIMENSION(nx,ny,ksta:kend) :: vout
    REAL(KIND=GP),INTENT (IN),DIMENSION(nx,ny,ksta:kend) :: vin
    INTEGER                                              :: i,j,k
!$omp parallel do if (kend-ksta.ge.nth) private (i,k)
    DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
      DO j = 1, ny
        DO i = 1, nx
          vout(i,j,k) = vin(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE R3toR3

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  METHOD     : Resize_Arrays
  !!  DESCRIPTION: Resize all arrays in the GPart class (including 
  !!               subclases, i.e. communicator, spline, and the
  !!               workspace), but not particle states.
  !!  ARGUMENTS  :
  !!    this    : 'this' class instance
  !!    new_size: new number of particles
  !!    onlyinc : if true, will only resize to increase array size
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ResizeArrays(this,new_size,onlyinc,exc)
!$  USE threads 
    IMPLICIT NONE
    CLASS(ParticleBase) ,INTENT(INOUT)             :: this
    INTEGER      ,INTENT(IN)                       :: new_size
    LOGICAL      ,INTENT(IN)                       :: onlyinc
    LOGICAL      ,INTENT(IN)          ,OPTIONAL    :: exc
    INTEGER                                        :: n
    LOGICAL                                        :: bret,assocptr(3)

    n = SIZE(this%id_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_IntArray(this%id_,new_size,.true.)
    END IF
    n = SIZE(this%tmpint_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_IntArray(this%tmpint_,new_size,.true.)
    END IF
    n = SIZE(this%ptmp0_,2)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank2(this%ptmp0_,new_size,.true.)
    END IF
    n = SIZE(this%gptmp0_,2)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank2(this%gptmp0_,new_size,.false.)
    END IF

    ! Resize workspace
    assocptr = .false.
    n = this%workspace_%get_nparts()
    IF ( ASSOCIATED(this%lvx_) ) THEN
      assocptr(1) = .true.
      call this%workspace_%free_pcomp_tmp(this%lvx_)
    END IF
    IF ( ASSOCIATED(this%lvy_) ) THEN
      assocptr(2) = .true.
      call this%workspace_%free_pcomp_tmp(this%lvy_)
    END IF
    IF ( ASSOCIATED(this%lvz_) ) THEN
      assocptr(3) = .true.
      call this%workspace_%free_pcomp_tmp(this%lvz_)
    END IF
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      call this%workspace_%resize_pcomp_arrays(new_size,.false.)
      call this%workspace_%set_nparts(new_size)
    END IF
    IF (assocptr(1)) call this%workspace_%get_pcomp_tmp(this%lvx_,bret)
    IF (assocptr(2)) call this%workspace_%get_pcomp_tmp(this%lvy_,bret)
    IF (assocptr(3)) call this%workspace_%get_pcomp_tmp(this%lvz_,bret)

    ! Resize VDB
    IF (this%iexchtype_.EQ.GPEXCHTYPE_VDB) THEN
      n = SIZE(this%vdb_)
      IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
        CALL Resize_ArrayRank2(this%vdb_,new_size,.false.)
      END IF
    ENDIF

    IF (PRESENT(exc)) THEN
      IF (exc) RETURN    ! Skip subclass resizing
    END IF

    CALL this%intop_ %ResizeArrays(new_size,onlyinc)
    CALL this%gpcomm_%ResizeArrays(new_size,onlyinc)
    RETURN 
  END SUBROUTINE ResizeArrays

end module particlebase_mod
