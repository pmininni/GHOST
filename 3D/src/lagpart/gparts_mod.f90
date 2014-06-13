!=================================================================
! GHOST GPart Lagrangian particles class
!
! 2013 D. Rosenberg
!      NCCS: ORNL
!
! 15 Jan 2013: Initial version
!=================================================================
MODULE class_GPart
      USE mpivars
      USE fprecision
      USE pdbtypes
      USE gtimer
      USE class_GPartComm
      USE class_GPSplineInt

      IMPLICIT NONE
      INCLUDE 'mpif.h' 

      INTEGER,PARAMETER,PUBLIC                       :: GPINIT_RANDLOC =0
      INTEGER,PARAMETER,PUBLIC                       :: GPINIT_USERLOC =1

      INTEGER,PARAMETER,PUBLIC                       :: GPINTRP_CSPLINE=0
      INTEGER,PARAMETER,PUBLIC                       :: GPINTRP_LAGINT =1

      INTEGER,PARAMETER,PUBLIC                       :: GPEXCHTYPE_NN  =0
      INTEGER,PARAMETER,PUBLIC                       :: GPEXCHTYPE_VDB =1

      INTEGER,PARAMETER,PUBLIC                       :: GPTIME_STEP    =1
      INTEGER,PARAMETER,PUBLIC                       :: GPTIME_COMM    =2
      INTEGER,PARAMETER,PUBLIC                       :: GPTIME_SPLINE  =3
      INTEGER,PARAMETER,PUBLIC                       :: GPTIME_TRANSP  =4
      INTEGER,PARAMETER,PUBLIC                       :: GPTIME_DATAEX  =5
      INTEGER,PARAMETER,PUBLIC                       :: GPTIME_INTERP  =6
      INTEGER,PARAMETER,PUBLIC                       :: GPTIME_PUPDATE =7

      INTEGER,PARAMETER,PRIVATE                      :: GPMAXTIMERS    =7  ! no. GPTIME parameters
      CHARACTER(len=8),PUBLIC                        :: lgext             ! string to hold time index
      CHARACTER(len=6),PUBLIC                        :: lgfmtext='(i8.8)' ! file time index format

      PRIVATE
      TYPE, PUBLIC :: GPart
        PRIVATE
        ! Member data:
        INTEGER, DIMENSION(MPI_STATUS_SIZE)          :: istatus_
        INTEGER                                      :: inittype_
        INTEGER                                      :: iinterp_
        INTEGER                                      :: iexchtype_
        INTEGER                                      :: iouttype_
        INTEGER                                      :: itimetype_
        TYPE(GPartComm)                              :: gpcomm_
        TYPE(GPSplineInt)                            :: intop_
        INTEGER                                      :: intorder_,itorder_,nd_(3),libnds_(3,2)
        INTEGER                                      :: myrank_,nprocs_
        INTEGER                                      :: htimers_(GPMAXTIMERS)
        INTEGER                                      :: ierr_,iseed_,istep_
        INTEGER                                      :: maxparts_,nparts_,nvdb_
        INTEGER                                      :: comm_
        INTEGER      , ALLOCATABLE, DIMENSION    (:) :: id_
!!      TYPE(GPDBrec), ALLOCATABLE, DIMENSION    (:) :: pdb_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: px_,py_,pz_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: lvx_,lvy_,lvz_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION  (:,:) :: ptmp0_,ptmp1_,vdb_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION  (:,:) :: gptmp0_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION    (:) :: ltmp0_,ltmp1_
        REAL(KIND=GP)                                :: lxbnds_(3,2),gext_(3)
        CHARACTER(len=1024)                          :: seedfile_,serr_,sfile_
      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: GPart_ctor
        FINAL            :: GPart_dtor
        PROCEDURE,PUBLIC :: Init            => GPart_Init
        PROCEDURE,PUBLIC :: Step            => GPart_StepRKK
        PROCEDURE,PUBLIC :: SetStep         => GPart_SetStepRKK
        PROCEDURE,PUBLIC :: SetLagVec       => GPart_SetLagVec
        PROCEDURE,PUBLIC :: FinalStep       => GPart_FinalizeRKK
        PROCEDURE,PUBLIC :: io_write_euler  => GPart_io_write_euler
        PROCEDURE,PUBLIC :: io_write_pdb    => GPart_io_write_pdb
        PROCEDURE,PUBLIC :: io_write_vec    => GPart_io_write_vec
        PROCEDURE,PUBLIC :: io_read         => GPart_io_read
        PROCEDURE,PUBLIC :: EulerToLag      => GPart_EulerToLag
        PROCEDURE,PUBLIC :: SetInitType     => GPart_SetInitType
        PROCEDURE,PUBLIC :: SetSeedFile     => GPart_SetSeedFile
        PROCEDURE,PUBLIC :: SetRandSeed     => GPart_SetRandSeed
        PROCEDURE,PUBLIC :: SetTimeOrder    => GPart_SetTimeOrder
        PROCEDURE,PUBLIC :: GetSeedFile     => GPart_GetSeedFile
        PROCEDURE,PUBLIC :: GetRandSeed     => GPart_GetRandSeed
        PROCEDURE,PUBLIC :: GetTimeOrder    => GPart_GetTimeOrder
        PROCEDURE,PUBLIC :: GetVDB          => GPart_GetVDB
        PROCEDURE,PUBLIC :: GetNParts       => GPart_GetNParts
        PROCEDURE,PUBLIC :: GetVel          => GPart_GetVel
        PROCEDURE,PUBLIC :: GetTime         => GPart_GetTime
        PROCEDURE,PUBLIC :: GetLoadBal      => GPart_GetLoadBal
!       PROCEDURE,PUBLIC :: GetPos
      END TYPE GPart

      PRIVATE :: GPart_Init            , GPart_StepRKK     
      PRIVATE :: GPart_SetStepRKK      , GPart_FinalizeRKK
      PRIVATE :: GPart_io_write_pdb    , GPart_io_read     
      PRIVATE :: GPart_io_write_euler  , GPart_io_write_vec
      PRIVATE :: GPart_InitRandSeed    , GPart_InitUserSeed 
      PRIVATE :: GPart_SetInitType     , GPart_SetSeedFile
      PRIVATE :: GPart_SetRandSeed     , GPart_Delete
      PRIVATE :: GPart_MakePeriodicP
      PRIVATE :: GPart_GetLocalWrk     , GPart_MakePeriodicExt
      PRIVATE :: GPart_GetLocalWrk_aux 
      PRIVATE :: GPart_ascii_write_pdb , GPart_binary_write_pdb
      PRIVATE :: GPart_ascii_write_lag , GPart_binary_write_lag
      PRIVATE :: GPart_ascii_read_pdb  , GPart_binary_read_pdb
      PRIVATE :: GPart_GetVDB          , GPart_GetVel
      PRIVATE :: GPart_GetTime         , GPart_GetLoadBal

! Methods:
  CONTAINS

  SUBROUTINE GPart_ctor(this,comm,mparts,inittype,iinterp,intorder,iexchtyp,iouttyp,csize,nstrip)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Main explicit constructor
!  ARGUMENTS:
!    this    : 'this' class instance
!    comm    : MPI communicator
!    mparts  : max no. particles allowed on grid. If inittype = GPINIT_USERLOC
!              then this is computed from the user file, and can
!              be any value.
!    inittype: GPINIT-typed quantity to give the type of particle
!              initialization
!    iinterp : GPINTRP-type quantity providing the interpolation
!              scheme
!    intorder: for variable-order (e.g., Lagrange) interpolation, 
!              the order (1, 2, 3...). Sets the number of 'ghost' zones
!              of data transferred between MPI tasks.
!    iexchtyp: format for exchanging particle data: GPEXCHTYPE_NN=nearest neighbor, 
!              suggesting that data is passed between neighboring
!              tasks only; GPEXCHTYPE_VDB form suggests that all tasks retain a copy
!              of the 'master' particle d.b., so that passing between
!              neighbors is not necessary. But the VDB form does require
!              a global reduction, and may be more expensive.
!    iouttup  : output type: 0==binary, 1=ASCII
!    csize    : cache size param for local transposes
!    nstrip   : 'strip-mining' size for local transposes
!-----------------------------------------------------------------
    USE grid
    USE mpivars
    USE commtypes
    USE random

    IMPLICIT NONE
    CLASS(GPart)     ,INTENT(INOUT)     :: this
    INTEGER          ,INTENT   (IN)     :: comm,mparts,csize,nstrip
    INTEGER                             :: disp(3),lens(3),types(3),szreal
    INTEGER          ,INTENT   (IN)     :: iexchtyp,iinterp,inittype,intorder,iouttyp
    INTEGER                             :: j,nc

    this%nparts_   = 0 
    this%nvdb_     = 0
    this%comm_     = comm
    this%maxparts_ = mparts
    this%nd_(1:3)  = n
    this%seedfile_ = 'gploc.dat'
    this%iinterp_  = 3          ! fixed for now
    this%inittype_ = inittype
    this%itorder_  = 2
    this%intorder_ = max(intorder,1)
    this%iseed_     = 1000
    this%istep_     = 0   
    this%iexchtype_=  iexchtyp
    this%iouttype_ =  iouttyp
    this%itimetype_=  GT_WTIME

    CALL prandom_seed(this%iseed_)
    IF ( this%intorder_ .NE. 3 ) THEN
      WRITE(*,*) 'GPart::ctor: Only 3rd order allowed for now' 
    ENDIF

    CALL MPI_COMM_SIZE(this%comm_,this%nprocs_,this%ierr_)
    CALL MPI_COMM_RANK(this%comm_,this%myrank_,this%ierr_)
    
    ! Iniitialze timers (get handles):
    DO j = 1, GPMAXTIMERS
!!    this%htimers_(j) = GTGetHandle()
      CALL GTStart(this%htimers_(j),this%itimetype_)
      IF ( this%htimers_(j).EQ.GTNULLHANDLE ) THEN
        WRITE(*,*) 'GPart_ctor: Not enough timers available'
        STOP
      ENDIF
    ENDDO

    CALL this%gpcomm_%GPartComm_ctor(GPCOMM_INTRFC_SF,this%maxparts_, &
         this%nd_,this%intorder_-1,this%comm_,this%htimers_(GPTIME_COMM))
    CALL this%gpcomm_%SetCacheParam(csize,nstrip)
    CALL this%gpcomm_%Init()

    DO j = 1,2
      this%libnds_(j,1) = 1 ; 
      this%libnds_(j,2) = n ; 
      this%lxbnds_(j,1) = 0.0_GP
      this%lxbnds_(j,2) = real(n-1,kind=GP)
    ENDDO
    this%libnds_(3,1) = ksta ; 
    this%libnds_(3,2) = kend ; 
    this%lxbnds_(3,1) = real(ksta-1,kind=GP) - 0.50_GP
    this%lxbnds_(3,2) = real(kend-1,kind=GP) + 0.50_GP !- 1.0_GP*epsilon(1.0_GP)
    IF ( this%myrank_ .EQ. 0 ) THEN
      this%lxbnds_(3,1) = -1.0_GP
    ENDIF
    IF ( this%myrank_ .EQ. this%nprocs_-1 ) THEN
      this%lxbnds_(3,2) = real(kend,kind=GP) 
    ENDIF

    DO j = 1,3
      this%gext_ (j) = real(this%nd_(j),kind=GP)
    ENDDO

    ! Instantiate interp operation. Remember that a valid timer 
    ! handle must be passed:
    CALL this%intop_%GPSplineInt_ctor(3,this%nd_,this%libnds_,this%lxbnds_,&
         this%maxparts_,this%gpcomm_,this%htimers_(GPTIME_DATAEX),&
         this%htimers_(GPTIME_TRANSP))

    ! Create part. d.b. structure type for I/O
    CALL MPI_TYPE_SIZE(GC_REAL,szreal,this%ierr_)

    ALLOCATE(this%id_      (this%maxparts_))
    ALLOCATE(this%px_      (this%maxparts_))
    ALLOCATE(this%py_      (this%maxparts_))
    ALLOCATE(this%pz_      (this%maxparts_))
    ALLOCATE(this%lvx_     (this%maxparts_))
    ALLOCATE(this%lvy_     (this%maxparts_))
    ALLOCATE(this%lvz_     (this%maxparts_))
    ALLOCATE(this%ptmp0_ (3,this%maxparts_))
    ALLOCATE(this%ptmp1_ (3,this%maxparts_))
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN
      ALLOCATE(this%gptmp0_ (3,this%maxparts_))
      ALLOCATE(this%vdb_(3,this%maxparts_))
    ENDIF
    ALLOCATE(this%ltmp0_ (this%maxparts_))
    ALLOCATE(this%ltmp1_ (this%maxparts_))


  END SUBROUTINE GPart_ctor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPart_dtor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Main explicit destructor
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    TYPE(GPart)   ,INTENT(INOUT)             :: this
    INTEGER                                  :: j


!!  CALL this%gpcomm_%GPartComm_dtor()

    IF ( ALLOCATED    (this%id_) ) DEALLOCATE   (this%id_)
    IF ( ALLOCATED    (this%px_) ) DEALLOCATE   (this%px_)
    IF ( ALLOCATED    (this%py_) ) DEALLOCATE   (this%py_)
    IF ( ALLOCATED    (this%pz_) ) DEALLOCATE   (this%pz_)
    IF ( ALLOCATED   (this%lvx_) ) DEALLOCATE  (this%lvx_)
    IF ( ALLOCATED   (this%lvy_) ) DEALLOCATE  (this%lvy_)
    IF ( ALLOCATED   (this%lvz_) ) DEALLOCATE  (this%lvz_)
    IF ( ALLOCATED (this%ptmp0_) ) DEALLOCATE(this%ptmp0_)
    IF ( ALLOCATED (this%gptmp0_) ) DEALLOCATE(this%gptmp0_)
    IF ( ALLOCATED (this%ptmp1_) ) DEALLOCATE(this%ptmp1_)
    IF ( ALLOCATED   (this%vdb_) ) DEALLOCATE  (this%vdb_)
    IF ( ALLOCATED (this%ltmp0_) ) DEALLOCATE(this%ltmp0_)
    IF ( ALLOCATED (this%ltmp1_) ) DEALLOCATE(this%ltmp1_)
    IF ( ALLOCATED   (this%lvy_) ) DEALLOCATE  (this%lvy_)
    IF ( ALLOCATED   (this%lvz_) ) DEALLOCATE  (this%lvz_)

    ! Destroy timers:
    DO j = 1, GPMAXTIMERS
      CALL GTFree(this%htimers_(j))
    ENDDO
  
  END SUBROUTINE GPart_dtor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPart_Init(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Init
!  DESCRIPTION: Initializes particle locations before integration.
!               Call after construction.
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart)   ,INTENT(INOUT)   :: this
    INTEGER                         :: j

    IF      ( this%inittype_ .EQ. GPINIT_RANDLOC ) THEN
      CALL GPart_InitRandSeed (this)   
    ELSE IF ( this%inittype_ .EQ. GPINIT_USERLOC ) THEN
      CALL GPart_InitUserSeed (this)
    ENDIF


  END SUBROUTINE GPart_Init
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_SetTimeOrder(this, iorder)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : SetTimeOrder
!  DESCRIPTION: Sets time step order
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iorder  : value of time step order
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)        :: this
    INTEGER      ,INTENT   (IN)        :: iorder

    this%itorder_ = iorder;
   
  END SUBROUTINE GPart_SetTimeOrder
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  FUNCTION GPart_GetTimeOrder(this) result(get_res)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetTimeOrder
!  DESCRIPTION: Gets time step order
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart),INTENT(INOUT)        :: this
    INTEGER                           :: get_res

    get_res = this%itorder_ 
   
  END FUNCTION GPart_GetTimeOrder
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_SetRandSeed(this, iseed)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : SetRandSeed
!  DESCRIPTION: Sets seed for random number generator
!  ARGUMENTS  :
!    this    : 'this' class instance
!    seed    : value of seed
!-----------------------------------------------------------------
    USE random

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)      :: this
    INTEGER      ,INTENT  (IN)       :: iseed

    this%iseed_     = iseed;
    CALL prandom_seed(this%iseed_)
    
   
  END SUBROUTINE GPart_SetRandSeed
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  FUNCTION GPart_GetRandSeed(this) result(get_res)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetRandSeed
!  DESCRIPTION: Gets random seed
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)       :: this
    INTEGER                           :: get_res

    get_res = this%iseed_
   
  END FUNCTION GPart_GetRandSeed
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_SetInitType(this, itype)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : SetInitType
!  DESCRIPTION: Sets particle seed initialization method
!  ARGUMENTS  :
!    this    : 'this' class instance
!    itype   : value of type GPINIT
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart),INTENT(INOUT)          :: this
    INTEGER     ,INTENT   (IN)          :: itype

    this%inittype_ = itype;
   
  END SUBROUTINE GPart_SetInitType
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  FUNCTION GPart_GetInitType(this) result(get_res)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetInitType
!  DESCRIPTION: Gets particle seed initialization method
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart),INTENT(INOUT)        :: this
    INTEGER                           :: get_res

    get_res = this%inittype_ 
   
  END FUNCTION GPart_GetInitType
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_SetSeedFile(this, sname)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : SetSeedFile
!  DESCRIPTION: Sets particle seed file name. Note: same file name
!               can be used for user-specified seeds, via InitUserSeed
!               or fixed-location seeds via InitFixedSeed.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    sname   : character file name
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart)    ,INTENT(INOUT)    :: this
    CHARACTER(len=*),INTENT   (IN)    :: sname

    this%seedfile_ = sname;
   
  END SUBROUTINE GPart_SetSeedFile
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  FUNCTION GPart_GetSeedFile(this) result(get_res)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetSeedFile
!  DESCRIPTION: Gets particle seed file name. 
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPart)     ,INTENT(INOUT)  :: this
    CHARACTER(1024)                  :: get_res

    get_res = this%seedfile_

  END FUNCTION GPart_GetSeedFile
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_InitRandSeed(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : InitRandSeed
!  DESCRIPTION: Initializes particle locations by dividing
!               maxparts evenly among MPI slabs, and randomly 
!               selecting positions within a slab.
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    USE random

    IMPLICIT NONE
    CLASS(GPart)  ,INTENT(INOUT)      :: this

    INTEGER                           :: ib,ie,j,iwrk1,iwrk2,nt
    REAL(KIND=GP)                     :: c,r,x1,x2
  
    ! Note: Box is [0,N-1]^3.
    ! All tasks write to _global_ grid, expecting a 
    ! VDBSynch afterwards, and a GetLocalWrk call to get the
    ! particles a task owns:
    iwrk1 = this%maxparts_/this%nprocs_
    iwrk2 = modulo(this%maxparts_,this%nprocs_)
    IF ( iwrk1.GT.0 ) THEN
      ib = this%myrank_*iwrk1+1+min(this%myrank_,iwrk2)
      ie = ib+iwrk1-1
      IF ( iwrk2.GT.this%myrank_ ) ie = ie + 1
    ELSE
      ib = this%myrank_+1
      ie = ib
      IF ( this%myrank_+1.GT.this%maxparts_ ) THEN
        ie = 0
        ib = 1
      ENDIF
    ENDIF
    ib = ib - 1
    ie = ie - 1
    this%nparts_ = ie - ib + 1
    DO j = 1, this%nparts_
       this%id_(j)    = ib + j - 1
       CALL prandom_number(r)
       x1 = real(this%libnds_(1,1)-1,kind=GP); x2 = real(this%libnds_(1,2)-1,kind=GP);
       c = r*(this%nd_(1)-1);
       this%px_(j) = min(max(c,x1),x2)
       CALL prandom_number(r)
       x1 = real(this%libnds_(2,1)-1,kind=GP); x2 = real(this%libnds_(2,2)-1,kind=GP);
       c = r*(this%nd_(1)-1);
       this%py_(j) = min(max(c,x1),x2)
       CALL prandom_number(r)
       x1 = real(this%libnds_(3,1)-1,kind=GP); x2 = real(this%libnds_(3,2)-1,kind=GP);
       this%pz_(j) = min(max(x1+r*(x2-x1),x1),x2)
    ENDDO
    CALL MPI_ALLREDUCE(this%nparts_,nt,1,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
    IF ( this%myrank_.eq.0 .AND. nt.NE.this%maxparts_ ) THEN
      WRITE(*,*) 'GPart_InitRandSeed: Inconsistent particle count: maxparts=', &
      this%maxparts_,' total created: ',nt
      STOP
    ENDIF
    CALL this%gpcomm_%VDBSynch(this%vdb_,this%maxparts_,this%id_, &
                          this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
    CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
                          this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
    CALL GPart_GetLocalWrk(this,this%id_,this%px_,this%py_,this%pz_,this%nparts_, &
                           this%vdb_,this%maxparts_)
    CALL GPART_ascii_write_pdb(this,1,'.','xlgInitRndSeed','000',0.0,this%vdb_)
    IF ( .NOT.GPart_PartNumConsistent(this,this%nparts_) ) THEN
      IF ( this%myrank_.eq.0 ) THEN
        WRITE(*,*) 'GPart_InitRandSeed: Invalid particle after GetLocalWrk call'
        STOP
      ENDIF
    ENDIF

  END SUBROUTINE GPart_InitRandSeed
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_InitUserSeed(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : InitUserSeed
!  DESCRIPTION: Initializes particle locations by letting user
!               specify particle number and locations in
!               file in member data seedfile_.
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    USE mpivars
    USE fprecision
    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)       :: this

    INTEGER                           :: navg,nl,nowned,nt
    INTEGER,ALLOCATABLE,DIMENSION(:)  :: iproc,ilproc
    REAL(KIND=GP)                     :: x,y,z


    ! Note: each record (line) consists of x y z real positions
    ! within [0,N-1]^3 box
    OPEN(UNIT=5,FILE=trim(this%seedfile_),STATUS='OLD',ACTION='READ',&
         IOSTAT=this%ierr_, IOMSG=this%serr_);
    IF ( this%ierr_ .NE. 0 ) THEN
      WRITE(*,*)'GPart::InitUserSeed: file:',this%seedfile_,' err: ', trim(this%serr_) 
      STOP
    ENDIF

    nt = 0  ! global part. record counter
    nl = 0  ! local particle counter
    DO WHILE ( this%ierr_.EQ.0 .AND. nt.LT.this%maxparts_ )
      READ(5,*,IOSTAT=this%ierr_) x, y, z
      IF ( this%ierr_ .NE. 0 ) EXIT
      IF ( z.GE.this%lxbnds_(3,1) .AND. z.LT.this%lxbnds_(3,2) .AND. &
           y.GE.this%lxbnds_(2,1) .AND. y.LT.this%lxbnds_(2,2) .AND. &
           x.GE.this%lxbnds_(1,1) .AND. x.LT.this%lxbnds_(1,2) ) THEN
        nl = nl + 1
        this%id_(nl) = nt
        this%px_(nl) = x
        this%py_(nl) = y
        this%pz_(nl) = z
      ENDIF
      nt = nt + 1
    ENDDO
    CLOSE(5)

    this%nparts_ = nl;
    CALL MPI_ALLREDUCE(nl,nt,1,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
    IF ( this%myrank_.eq.0 .AND. nt.NE.this%maxparts_ ) THEN
      WRITE(*,*) 'GPart_InitUserSeed: Inconsistent particle count: maxparts=', &
      this%maxparts_,' total read: ',nt
      STOP
    ENDIF
    CALL this%gpcomm_%VDBSynch(this%vdb_,this%maxparts_,this%id_, &
                          this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)


  END SUBROUTINE GPart_InitUserSeed
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_io_write_pdb(this, iunit, dir, spref, nmb, time)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : io_write_pdb
!  DESCRIPTION: Does write of Lagrangian position d.b. to file. 
!               Position of the particle structure in file is the
!               particle's id. Main entry point for both ASCII and
!               binary writes.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iunit   : unit number
!    dir     : output directory
!    spref   : filename prefix
!    nmd     : time index
!    time    : real time
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)       :: this
    REAL(KIND=GP),INTENT   (IN)       :: time
    REAL(KIND=GP)                     :: prec(3)
    INTEGER,INTENT(IN)                :: iunit
    INTEGER                           :: fh,j,nt
    INTEGER(kind=MPI_OFFSET_KIND)     :: offset
    CHARACTER(len=*),INTENT(IN)       :: dir
    CHARACTER(len=*),INTENT(IN)       :: nmb
    CHARACTER(len=*),INTENT(IN)       :: spref
    TYPE(GPDBrec)                     :: pst

    ! Do a sanity check:
!!  CALL MPI_ALLREDUCE(this%nparts_,nt,1,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
!!  IF ( nt .NE. this%maxparts_ ) THEN
!!    WRITE(*,*) this%myrank_, ': GPart_io_write_pdb: particle inconsistency: no. required=',&
!!               this%maxparts_,' no. found=',nt
!!    STOP
!!  ENDIF

    IF ( this%iexchtype_.EQ.GPEXCHTYPE_NN ) THEN
      CALL this%gpcomm_%VDBSynch(this%ptmp0_,this%maxparts_,this%id_, &
           this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
    ELSE
      Do j = 1, this%maxparts_
        this%ptmp0_(1:3,j) = this%vdb_(1:3,j)
      ENDDO
    ENDIF

    IF ( this%iouttype_ .EQ. 0 ) THEN
      CALL GPart_binary_write_pdb(this,iunit,dir,spref,nmb,time,this%ptmp0_)
    ELSE
      CALL GPart_ascii_write_pdb (this,iunit,dir,spref,nmb,time,this%ptmp0_)
    ENDIF

  END SUBROUTINE GPart_io_write_pdb


  SUBROUTINE GPart_io_write_vec(this, iunit, dir, spref, nmb, time)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : io_write_vec
!  DESCRIPTION: Does write of Lagrangian vector that is
!               currently stored. This vector may not be the
!               advecting velocities if a call to SetLagVec
!               is made with a different vector.
!               
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iunit   : unit number
!    dir     : output directory
!    spref   : filename prefix
!    nmd     : time index
!    time    : real time
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)       :: this
    REAL(KIND=GP),INTENT   (IN)       :: time
    INTEGER,INTENT(IN)                :: iunit
    CHARACTER(len=*),INTENT(IN)       :: dir
    CHARACTER(len=*),INTENT(IN)       :: nmb
    CHARACTER(len=*),INTENT(IN)       :: spref
    TYPE(GPDBrec)                     :: pst

    IF ( .NOT.GPart_PartNumConsistent(this,this%nparts_) ) THEN
      IF ( this%myrank_.eq.0 ) THEN
        WRITE(*,*) 'GPart_io_write_vec: Inconsistent particle count'
        STOP
      ENDIF
    ENDIF

!   CALL this%gpcomm_%VDBSynch(this%ptmp0_,this%maxparts_,this%id_, &
!        this%lvx_,this%lvy_,this%lvz_,this%nparts_,this%ptmp1_)

    IF ( this%iouttype_ .EQ. 0 ) THEN
      CALL GPart_binary_write_lag(this,iunit,dir,spref,nmb,time,this%lvx_,this%lvy_,this%lvz_)
    ELSE
      !CALL GPart_ascii_write_lag(this,iunit,dir,spref,nmb,time,this%lvx_,this%lvy_,this%lvz_)
      CALL this%gpcomm_%VDBSynch(this%ptmp0_,this%maxparts_,this%id_, &
           this%lvx_,this%lvy_,this%lvz_,this%nparts_,this%ptmp1_)
      CALL GPart_ascii_write_lag(this,iunit,dir,spref,nmb,time, &
                                 this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:))
    ENDIF

  END SUBROUTINE GPart_io_write_vec
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_binary_write_pdb(this, iunit, dir, spref, nmb, time, pdb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : binary write
!  DESCRIPTION: Does binary write of Lagrangian position d.b. to file. 
!               Position of the particle structure in file is the
!               particle's id.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iunit   : unit number
!    dir     : output directory
!    spref   : filename prefix
!    nmd     : time index
!    time    : real time
!    pdb     : particle d.b. in (3,maxparts) array
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)       :: this
    REAL(KIND=GP),INTENT   (IN)       :: time
    REAL(KIND=GP),INTENT   (IN)       :: pdb(3,this%maxparts_)
    REAL(KIND=GP)                     :: prec(3)
    INTEGER,INTENT(IN)                :: iunit
    INTEGER                           :: fh,nt,szreal
    INTEGER(kind=MPI_OFFSET_KIND)     :: offset
    CHARACTER(len=*),INTENT(IN)       :: dir
    CHARACTER(len=*),INTENT(IN)       :: nmb
    CHARACTER(len=*),INTENT(IN)       :: spref
    TYPE(GPDBrec)                     :: pst

    INTEGER                           :: j,gc,lc

    ! Must write part. data to correct position in file:
    CALL MPI_TYPE_SIZE(GC_REAL    ,szreal,this%ierr_)
    CALL MPI_FILE_OPEN(this%comm_,trim(dir) // '/' // trim(spref) // &
         '.' // nmb // '.lag',MPI_MODE_CREATE+MPI_MODE_WRONLY, &
          MPI_INFO_NULL,fh,this%ierr_)
    prec(1) = real(this%maxparts_,kind=GP)
    offset = 0
    CALL MPI_FILE_WRITE_AT_ALL(fh,offset,prec(1),1,GC_REAL,this%istatus_,this%ierr_)
    offset = szreal
    CALL MPI_FILE_WRITE_AT_ALL(fh,offset,time   ,1,GC_REAL,this%istatus_,this%ierr_)
    gc = 0
    DO j = 1, this%nparts_
      offset  = (this%id_(j)*3+2)*szreal
      prec(1) = this%px_(j)
      prec(2) = this%py_(j)
      prec(3) = this%pz_(j)
      CALL MPI_FILE_WRITE_AT(fh,offset,prec,3,GC_REAL,this%istatus_,this%ierr_)
      CALL MPI_GET_COUNT(this%istatus_,GC_REAL,lc,this%ierr_)
      gc = gc+lc
    ENDDO
    CALL MPI_FILE_CLOSE(fh,this%ierr_)

    IF ( gc .NE. this%nparts_*3 ) THEN
      WRITE(*,*)this%myrank_, &
        ': GPart_binary_write_pdb: insufficient amount of data written; no. required=',&
        this%nparts_*3,' no. written=',gc
      STOP
    ENDIF

  END SUBROUTINE GPart_binary_write_pdb
!-----------------------------------------------------------------
!-----------------------------------------------------------------

 SUBROUTINE GPart_ascii_write_pdb(this, iunit, dir, spref, nmb, time,pdb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : ascii_write_pdb
!  DESCRIPTION: Does ASCII write of Lagrangian position d.b. to file.
!               The local MPI tasks write to a file with prefix
!               spref, in the following format:
!                     dir/spref.TTT.PPP.txt
!               where TTT is the time index, given by nmb, and
!               PPP is the MPI rank.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iunit   : unit number
!    dir     : output directory
!    spref   : filename prefix
!    nmd     : time index
!    time    : real time
!    pdb     : part. d.b. in (4,maxparts) array
!
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)       :: this
    REAL(KIND=GP),INTENT   (IN)       :: time
    REAL(KIND=GP),INTENT   (IN)       :: pdb(3,this%maxparts_)
    INTEGER      ,INTENT   (IN)       :: iunit
    INTEGER                           :: j
    CHARACTER(len=*),INTENT(IN)       :: dir
    CHARACTER(len=*),INTENT(IN)       :: nmb
    CHARACTER(len=*),INTENT(IN)       :: spref
    CHARACTER(len=3)                  :: sind

      ! Write global VDB, with time header, indexed only
      ! by time index: dir/spref.TTT.txt:
     IF ( this%myrank_.EQ.0 ) THEN
       OPEN(iunit,file=trim(dir)// '/' // trim(spref) // '.' // &
                 nmb //  '.txt')
       WRITE(iunit,*) this%maxparts_
       WRITE(iunit,*) time
       DO j = 1, this%maxparts_
         WRITE(iunit,600) pdb(1,j),pdb(2,j),pdb(3,j)
  600    FORMAT(3(E23.15,1X))
       ENDDO
      CLOSE(iunit)
    ENDIF

  END SUBROUTINE GPart_ascii_write_pdb
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_io_write_euler(this, iunit, dir, spref, nmb, &
             time, evar, doupdate, tmp1, tmp2)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : io_write_euler
!  DESCRIPTION: Converts specified Eulerian real-space variable to
!               a Lagrangian quantity by interpolating to particle positions;
!               does write of Lagrangian variable to file. 
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iunit   : unit number
!    dir     : output directory
!    fname   : filename prefix
!    nmb     : time index
!    time    : real time
!    evar    : Eulerian data from which to compute Lagrangian 
!              quantity: theta(y) = theta(x(y),t). Interpolation
!              of evar is done internally before write. Note that
!              data in evar is lost on exit.
!    doupdate: if true, do interp point update in interpolator; else don't
!    tmp1/2  : temp arrays of size of evar. Required for interpolation
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)                         :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(n,n,ksta:kend):: evar(n,n,ksta:kend)
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(n,n,ksta:kend):: tmp1,tmp2
    REAL(KIND=GP),INTENT   (IN)                         :: time
    INTEGER      ,INTENT   (IN)                         :: iunit
    INTEGER                                             :: fh,offset,nt,szint,szreal
    INTEGER                                             :: j
    LOGICAL      ,INTENT   (IN)                         :: doupdate
    CHARACTER(len=100), INTENT(IN)                      :: dir
    CHARACTER(len=*)  , INTENT(IN)                      :: nmb
    CHARACTER(len=*)  , INTENT(IN)                      :: spref
    CHARACTER(len=1024)                                 :: sfile


    CALL GPart_EulerToLag(this,this%ltmp1_,this%nparts_,evar,doupdate,tmp1,tmp2)

    IF ( this%iouttype_ .EQ. 0 ) THEN
      CALL GPart_binary_write_lag(this,iunit,dir,spref,nmb,time,this%ltmp1_)
    ELSE
      ! First, must synch up the field data since only task 0 writes:
      CALL this%gpcomm_%LagSynch(this%ltmp0_,this%maxparts_,this%id_, &
                                 this%ltmp1_,this%nparts_,this%ptmp0_)
      CALL GPart_ascii_write_lag (this,iunit,dir,spref,nmb,time,this%ltmp0_)
    ENDIF

    
  END SUBROUTINE GPart_io_write_euler
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_binary_write_lag(this, iunit, dir, spref, nmb, time, &
             fld0, fld1, fld2, fld3, fld4, fld5, fld6, fld7, fld8)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPart_binary_write_lag
!  DESCRIPTION: Does binary write of Lagrangian field to file. 
!               Position of the particle structure in file is the
!               particle's id. This method allows for up to 9
!               Lagranian variables to be outputted. At least one
!               variable _must_ be present (fld0). Do not use keywords
!               to specify optional arguments. 
!
!               Note that this call will have all MPI tasks write
!               their data collectively, so no 'synching' of data 
!               is required on unput.
!
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iunit   : unit number
!    dir     : output directory
!    spref   : filename prefix
!    nmd     : time index
!    time    : real time
!    fld0-8  : Lagrangian field
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)          :: this
    REAL(KIND=GP),INTENT   (IN)          :: time
    REAL(KIND=GP),INTENT   (IN)          :: fld0(this%maxparts_)
    REAL(KIND=GP),INTENT   (IN),OPTIONAL,&
    DIMENSION(this%maxparts_)            :: fld1,fld2,fld3,fld4,fld5,fld6,fld7,fld8
    REAL(KIND=GP)                        :: vout(9)
    INTEGER,INTENT(IN)                   :: iunit
    INTEGER                              :: fh,nt,nv,szreal
    INTEGER(kind=MPI_OFFSET_KIND)        :: offset
    CHARACTER(len=*),INTENT(IN)          :: dir
    CHARACTER(len=*),INTENT(IN)          :: nmb
    CHARACTER(len=*),INTENT(IN)          :: spref
    TYPE(GPDBrec)                        :: pst

    INTEGER                              :: j,gc,lc

    nv = 1 
    IF ( present(fld1) ) nv=nv+1
    IF ( present(fld2) ) nv=nv+1
    IF ( present(fld3) ) nv=nv+1
    IF ( present(fld4) ) nv=nv+1
    IF ( present(fld5) ) nv=nv+1
    IF ( present(fld6) ) nv=nv+1
    IF ( present(fld7) ) nv=nv+1
    IF ( present(fld8) ) nv=nv+1

    ! Must write part. data to correct position in file:
    CALL MPI_TYPE_SIZE(GC_REAL    ,szreal,this%ierr_)
    CALL MPI_FILE_OPEN(this%comm_,trim(dir) // '/' // trim(spref) // &
         '.' // nmb // '.lag',MPI_MODE_CREATE+MPI_MODE_WRONLY, &
          MPI_INFO_NULL,fh,this%ierr_)
    offset = 0
    CALL MPI_FILE_WRITE_AT_ALL(fh,offset,real(this%maxparts_,kind=GP),1,GC_REAL,this%istatus_,this%ierr_)
    offset = szreal
    CALL MPI_FILE_WRITE_AT_ALL(fh,offset,time   ,1,GC_REAL,this%istatus_,this%ierr_)
    gc = 0
    DO j = 1, this%nparts_
      offset  = (nv*this%id_(j)+2)*szreal
      vout(1) = fld0(j)
      IF ( present(fld1) ) vout(2) = fld1(j)
      IF ( present(fld2) ) vout(3) = fld2(j)
      IF ( present(fld3) ) vout(4) = fld3(j)
      IF ( present(fld4) ) vout(5) = fld4(j)
      IF ( present(fld5) ) vout(6) = fld5(j)
      IF ( present(fld6) ) vout(7) = fld6(j)
      IF ( present(fld7) ) vout(8) = fld7(j)
      IF ( present(fld8) ) vout(9) = fld8(j)
      CALL MPI_FILE_WRITE_AT(fh,offset,vout,nv,GC_REAL,this%istatus_,this%ierr_)
      CALL MPI_GET_COUNT(this%istatus_,GC_REAL,lc,this%ierr_)
      gc = gc+lc
    ENDDO
    CALL MPI_FILE_CLOSE(fh,this%ierr_)

    IF ( gc .NE. this%nparts_*nv ) THEN
      WRITE(*,*)this%myrank_, &
        ': GPart_binary_write_lag: insufficient amount of data written; no. required=',&
        this%nparts_*nv,' no. written=',gc
      STOP
    ENDIF

  END SUBROUTINE GPart_binary_write_lag
!-----------------------------------------------------------------
!-----------------------------------------------------------------


 SUBROUTINE GPart_ascii_write_lag(this, iunit, dir, spref, nmb, time, &
            fld0, fld1, fld2, fld3, fld4, fld5, fld6, fld7, fld8)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPart_ascii_write_lag
!  DESCRIPTION: Does ASCII write of Lagrangian fld to file.
!               The local MPI tasks write to a file with prefix
!               spref, in the following format:
!                     dir/spref.TTT.PPP.txt
!               where TTT is the time index, given by nmb, and
!               PPP is the MPI rank.  This method allows for up to 9
!               Lagranian variables to be outputted. At least one
!               variable _must_ be present (fld0). Do not use keywords
!               to specify optional arguments.
!
!               Note that this call will have ony MPI task 0 write
!               data, so data should be 'synched' before calling. 
!               This is not done here.
!
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iunit   : unit number
!    dir     : output directory
!    spref   : filename prefix
!    nmd     : time index
!    time    : real time
!    fld0-8  : Lagrangian fields
!
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)       :: this
    REAL(KIND=GP),INTENT   (IN)       :: time
    REAL(KIND=GP),INTENT   (IN)       :: fld0(this%maxparts_)
    REAL(KIND=GP),INTENT   (IN),OPTIONAL,&
    DIMENSION(this%maxparts_)         :: fld1,fld2,fld3,fld4,fld5,fld6,fld7,fld8
    REAL(KIND=GP)                     :: vout(9)
    INTEGER      ,INTENT   (IN)       :: iunit
    INTEGER                           :: j,nv
    CHARACTER(len=*),INTENT(IN)       :: dir
    CHARACTER(len=*),INTENT(IN)       :: nmb
    CHARACTER(len=*),INTENT(IN)       :: spref
    CHARACTER(len=3)                  :: sind

    nv = 1 
    IF ( present(fld1) ) nv=nv+1
    IF ( present(fld2) ) nv=nv+1
    IF ( present(fld3) ) nv=nv+1
    IF ( present(fld4) ) nv=nv+1
    IF ( present(fld5) ) nv=nv+1
    IF ( present(fld6) ) nv=nv+1
    IF ( present(fld7) ) nv=nv+1
    IF ( present(fld8) ) nv=nv+1

    ! Write global VDB, with time header, indexed only
    ! by time index: dir/spref.TTT.txt:
    IF ( this%myrank_.EQ.0 ) THEN
      OPEN(iunit,file=trim(dir)// '/' // trim(spref) // '.' // &
            nmb //  '.txt')
      WRITE(iunit,*) this%maxparts_
      WRITE(iunit,*) time
      DO j = 1, this%maxparts_
        vout(1) = fld0(j)
        IF ( present(fld1) ) vout(2) = fld1(j)
        IF ( present(fld2) ) vout(3) = fld2(j)
        IF ( present(fld3) ) vout(4) = fld3(j)
        IF ( present(fld4) ) vout(5) = fld4(j)
        IF ( present(fld5) ) vout(6) = fld5(j)
        IF ( present(fld6) ) vout(7) = fld6(j)
        IF ( present(fld7) ) vout(8) = fld7(j)
        IF ( present(fld8) ) vout(9) = fld8(j)
        WRITE(iunit,600) vout(1:nv)
  600   FORMAT(9(E23.15,1X))
      ENDDO
     CLOSE(iunit)
   ENDIF

  END SUBROUTINE GPart_ascii_write_lag
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_io_read(this, iunit, dir, spref, nmb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : io_read
!  DESCRIPTION: Does read of Lagrangian particle data from file. 
!               This is the main entry point for both binary and
!               ASCII reads.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iunit   : unit number
!    dir     : input directory
!    spref   : filename prefix
!    nmb     : time index
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart)    ,INTENT(INOUT)            :: this
!!  REAL(KIND=GP),INTENT(OUT),DIMENSION(:)    :: lvar
    REAL(KIND=GP)                             :: rvar,time
    INTEGER,INTENT(IN)                        :: iunit
    INTEGER                                   :: fh,j,ng
    INTEGER(kind=MPI_OFFSET_KIND)             :: offset
    CHARACTER(len=*),INTENT   (IN)            :: dir
    CHARACTER(len=*),INTENT   (IN)            :: nmb
    CHARACTER(len=*),INTENT   (IN)            :: spref

    IF ( this%iouttype_ .EQ. 0 ) THEN
      CALL GPart_binary_read_pdb(this,iunit,dir,spref,nmb,time,this%ptmp0_)
    ELSE
      CALL GPart_ascii_read_pdb (this,iunit,dir,spref,nmb,time,this%ptmp0_)
    ENDIF

    CALL GPart_GetLocalWrk(this,this%id_,this%px_,this%py_,this%pz_, &
                           this%nparts_,this%ptmp0_,this%maxparts_)

    CALL MPI_ALLREDUCE(this%nparts_,ng,1,MPI_INTEGER,   &
                       MPI_SUM,this%comm_,this%ierr_)
    IF ( this%myrank_.EQ.0 .AND. ng.NE.this%maxparts_ ) THEN
      WRITE(*,*)'GPart_io_read: inconsistent d.b.: expected: ', &
                 this%maxparts_, '; found: ',ng
      STOP
    ENDIF


    ! If there is a global VDB for data 'exchanges', create it here:
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN
      DO j = 1, this%maxparts_
        this%vdb_(1:3,j) = this%ptmp0_(1:3,j)
      ENDDO
    ENDIF


  END SUBROUTINE GPart_io_read


  SUBROUTINE GPart_binary_read_pdb(this,iunit,dir,spref,nmb,time,pdb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : binary_read_pdb
!  DESCRIPTION: Does read of binary Lagrangian particle data from file. 
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iunit   : unit number
!    dir     : input directory
!    spref   : filename prefix
!    nmb     : time index
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)               :: this
    REAL(KIND=GP)                             :: rvar,time
    REAL(KIND=GP),INTENT(INOUT)               :: pdb(3,this%maxparts_)
    INTEGER,INTENT(IN)                        :: iunit
    INTEGER                                   :: fh,j,lc,szreal
    INTEGER(kind=MPI_OFFSET_KIND)             :: offset
    CHARACTER(len=*),INTENT   (IN)            :: dir
    CHARACTER(len=*),INTENT   (IN)            :: nmb
    CHARACTER(len=*),INTENT   (IN)            :: spref

    CALL MPI_TYPE_SIZE(GC_REAL,szreal,this%ierr_)
    CALL MPI_FILE_OPEN(this%comm_,trim(dir) // '/' // trim(spref) // &
         '.' // nmb // '.lag',MPI_MODE_RDONLY,MPI_INFO_NULL,fh,this%ierr_)
  
    ! Must read part. data from correct spot in file:
    offset = 0
    CALL MPI_FILE_READ_AT_ALL(fh,offset,rvar,1,GC_REAL,this%istatus_,this%ierr_)    !  no.parts
    IF ( int(rvar).NE.this%maxparts_ ) THEN
      WRITE(*,*) 'GPart_io_read: Attempt to read incorrect number of particles: required:',&
                  this%maxparts_,' no attempted: ',int(rvar)
      STOP
    ENDIF
    offset = szreal
    CALL MPI_FILE_READ_AT_ALL(fh, offset,rvar,1,GC_REAL,this%istatus_,this%ierr_) ! time
    offset = 2*szreal
    CALL MPI_FILE_READ_AT_ALL(fh,offset,pdb,3*this%maxparts_,GC_REAL,this%istatus_,this%ierr_)
    CALL MPI_FILE_CLOSE(fh,this%ierr_)

  END SUBROUTINE GPart_binary_read_pdb
!-----------------------------------------------------------------
!-----------------------------------------------------------------

 SUBROUTINE GPart_ascii_read_pdb(this,iunit,dir,spref,nmb,time,pdb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : ascii_read_pdb
!  DESCRIPTION: Does ASCII read of Lagrangian position d.b. from file.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iunit   : unit number
!    dir     : output directory
!    spref   : filename prefix
!    nmd     : time index
!    time    : real time
!    pdb     : part. d.b. in (4,maxparts) array
!
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)       :: this
    REAL(KIND=GP),INTENT(INOUT)       :: time
    REAL(KIND=GP),INTENT(INOUT)       :: pdb(3,this%maxparts_)
    INTEGER      ,INTENT   (IN)       :: iunit
    INTEGER                           :: j,kstat,nt
    CHARACTER(len=*),INTENT(IN)       :: dir
    CHARACTER(len=*),INTENT(IN)       :: nmb
    CHARACTER(len=*),INTENT(IN)       :: spref
    CHARACTER(len=3)                  :: sind

    ! Read global VDB, with time header, indexed only
    ! by time index: dir/spref.TTT.txt:
    IF ( this%myrank_.EQ.0 ) THEN
      OPEN(iunit,file=trim(dir)// '/' // trim(spref) // '.' // &
                nmb //  '.txt',status='old',form='formatted',iostat=kstat)
      IF ( kstat.NE.0 ) THEN
        WRITE(*,*)'GPart_ascii_read_pdb: could not open file for reading: ',&
        trim(dir)// '/' // trim(spref) // '.' // nmb //  '.txt'
        STOP
      ENDIF
      READ(iunit,*,iostat=kstat) nt
      READ(iunit,*,iostat=kstat) time
      IF ( nt.NE.this%maxparts_ ) THEN
        WRITE(*,*)this%myrank_, &
          ': GPart_ascii_read_pdb: particle inconsistency: no. required=',&
          this%maxparts_,' no. found=',nt, &
          ' file=',trim(dir)// '/' // trim(spref) // '.' // nmb //  '.txt'
        STOP
      ENDIF
      DO j = 1, this%maxparts_
        READ(iunit,*,iostat=kstat) pdb(1,j),pdb(2,j),pdb(3,j)
  600   FORMAT(3(E23.15,1X))
      ENDDO
      CLOSE(iunit)
    ENDIF
    CALL MPI_BCAST(pdb,3*this%maxparts_,GC_REAL,0,this%comm_,this%ierr_)
    IF ( this%ierr_.NE.MPI_SUCCESS ) THEN
        WRITE(*,*)this%myrank_, ': GPart_ascii_read_pdb: Broadcast failed: file=',&
        trim(dir)// '/' // trim(spref) // '.' // nmb //  '.txt'
    ENDIF

  END SUBROUTINE GPart_ascii_read_pdb
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_SetStepRKK(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     :SetStepRKK
!  DESCRIPTION: Initializes an explicit integration timstep.
!               Must be called at the start of the RK stage execution,.

!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart)   ,INTENT(INOUT)          :: this

    INTEGER                                :: j

    ! Initialize solution, u: 
    ! u* <-- u: 
 
    ! Cycle over JST loop to update state:
!$omp parallel do
    DO j = 1, this%nparts_
       this%ptmp0_(1,j) = this%px_(j)  ! u_0
       this%ptmp0_(2,j) = this%py_(j)  ! u_0
       this%ptmp0_(3,j) = this%pz_(j)  ! u_0
    ENDDO 

  END SUBROUTINE GPart_SetStepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_FinalizeRKK(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : FinalizeRKK
!  DESCRIPTION: Called at the end of all RK-like stages to
!               complete Lagrangian particle update.

!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)                 :: this
    INTEGER                                     :: j,ng

    ! u(t+dt) = u*: done already

    ! IF using nearest-neighbor interface, do particle exchange 
    ! between nearest-neighbor tasks BEFORE z-PERIODIZING particle coordinates:
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_NN ) THEN
      CALL GTStart(this%htimers_(GPTIME_COMM))
      CALL this%gpcomm_%PartExchangeV(this%id_,this%px_,this%py_,this%pz_, &
           this%nparts_,this%lxbnds_(3,1),this%lxbnds_(3,2))
      CALL GTAcc(this%htimers_(GPTIME_COMM))
    ENDIF

    ! Enforce periodicity in x, y, & z:
    CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,7)

    ! If using VDB interface, do synch-up, and get local work:
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN

      IF ( .NOT.GPart_PartNumConsistent(this,this%nparts_) ) THEN
        IF ( this%myrank_.eq.0 ) THEN
          WRITE(*,*) 'GPart_FinalizeRKK: Inconsistent particle count'
        ENDIF
      ENDIF
      ! Synch up VDB, if necessary:
      CALL GTStart(this%htimers_(GPTIME_COMM))
      CALL this%gpcomm_%VDBSynch(this%vdb_,this%maxparts_,this%id_, &
           this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
      CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
                     this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:),&
                     this%nparts_,this%ptmp1_)
      CALL GTAcc(this%htimers_(GPTIME_COMM))

      ! If using VDB, get local particles to work on:
      ! GPart_GetLocalWrk_aux also synchronizes auxiliary RK arrays, 
      ! and is needed if the call is done in the middle of a RK step.
      CALL GPart_GetLocalWrk_aux(this,this%id_,this%px_,this%py_,this%pz_,&
                       this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:),&
                       this%nparts_,this%vdb_,this%gptmp0_,this%maxparts_)
      CALL MPI_ALLREDUCE(this%nparts_,ng,1,MPI_INTEGER,   &
                         MPI_SUM,this%comm_,this%ierr_)

      IF ( this%myrank_.EQ.0 .AND. ng.NE.this%maxparts_) THEN
        WRITE(*,*)'GPart_FinalizeRKK: inconsistent d.b.: expected: ', &
                 this%maxparts_, '; found: ',ng
        CALL GPART_ascii_write_pdb(this,1,'.','xlgerr','000',0.0,this%vdb_)
        STOP
      ENDIF

    ENDIF
    
    RETURN

  END SUBROUTINE GPart_FinalizeRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_StepRKK(this, vx, vy, vz, dt, xk, tmp1, tmp2)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : StepRKK
!  DESCRIPTION: Carries out one stage of explicit RK-like time
!               inegration step.  Intended for explicit step within 
!               an outer stepper method of the form:
!
!               X = X_0 + dt * V(X(t),t) * xk,
!       
!               Note that the vx, vy, vz, will be overwritten here.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    vz,vy,vz: compoments of velocity field, in real space, partially
!              updated, possibly. These will be overwritten!
!    dt      : integration timestep
!    xk      : multiplicative RK time stage factor
!    tmp1(2) : temp arrays the same size as vx, vy, vz
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)                          :: this
    INTEGER                                              :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(n,n,ksta:kend) :: vx,vy,vz,tmp1,tmp2
    REAL(KIND=GP),INTENT   (IN)                          :: dt,xk
    REAL(KIND=GP)                                        :: dtfact
    REAL(KIND=GP),ALLOCATABLE  ,DIMENSION            (:) :: lid,gid

    dtfact = dt*xk*real(n,kind=GP)/(8.0_GP*atan(1.0_GP))
    CALL GTStart(this%htimers_(GPTIME_STEP))

    ! Find F(u*):
    ! ... x:
    CALL GPart_EulerToLag(this,this%lvx_,this%nparts_,vx,.true.,tmp1,tmp2)
    ! ux* <-- ux + dt * F(U*)*xk:
    DO j = 1, this%nparts_
      this%px_(j) = this%ptmp0_(1,j) + dtfact*this%lvx_(j)
    ENDDO
    !
    ! ... y:
    ! Exchange bdy data for velocities, so that we
    ! can perform local interpolations:
    CALL GPart_EulerToLag(this,this%lvy_,this%nparts_,vy,.false.,tmp1,tmp2)
    ! uy* <-- uy + dt * F(U*)*xk:
    DO j = 1, this%nparts_
      this%py_(j) = this%ptmp0_(2,j) + dtfact*this%lvy_(j)
    ENDDO

    ! ... z:
    ! Exchange bdy data for velocities, so that we
    ! can perform local interpolations:
    CALL GPart_EulerToLag(this,this%lvz_,this%nparts_,vz,.false.,tmp1,tmp2)
    ! uz* <-- uz + dt * F(U*)*xk:
    DO j = 1, this%nparts_
      this%pz_(j) = this%ptmp0_(3,j) + dtfact*this%lvz_(j)
    ENDDO

    ! Enforce periodicity in x-y only:
    CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,3)

    CALL GTAcc(this%htimers_(GPTIME_STEP))
!   ALLOCATE  (lid(this%maxparts_))
!   ALLOCATE  (gid(this%maxparts_))
!   lid = 0
!   gid = 0
!   do j=1,this%maxparts_
!     lid(this%id_(j)+1) = 1
!   enddo
!   call mpi_allreduce(lid,gid,this%maxparts_,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
!   if ( this%myrank_.eq. 0 ) then
!     do j=1,this%maxparts_
!       if ( gid(j) .gt.1 ) then
!         write(*,*)this%myrank_,': StepRKK: multiplicity > 1: id=',j-1,' p=', &
!                   this%vdb_(1,j),this%vdb_(2,j),this%vdb_(3,j)
!       else if ( gid(j).eq.0 ) then
!         write(*,*)this%myrank_,': StepRKK: particle missing: id=',j-1
!       endif    
!     enddo
!   endif
!
!   DEALLOCATE(lid,gid)

  END SUBROUTINE GPart_StepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------



  SUBROUTINE GPart_SetLagVec(this, vx, vy, vz, doupdate, tmp1, tmp2)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : SetLagVec
!  DESCRIPTION: Sets internal vector at Lagrangian particle positions.
!               This interface allows us to store the vectory in
!               the internal velocity vector locations, for e.g., output.
!               On entry, the interpolator is updated when this interface 
!               is called, only if doupdate is set to true. This update is 
!               done only for the first component, as it's not needed after this. 
!
!  ARGUMENTS  :
!    this    : 'this' class instance
!    vz,vy,vz: compoments of vector field, in real space, partially
!              updated, possibly. These will be overwritten!
!    doupdate: if true, do interp point update in interpolator; else don't
!    tmp1(2) : temp arrays the same size as vx, vy, vz
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)                          :: this
    LOGICAL      ,INTENT   (IN)                          :: doupdate
    INTEGER                                              :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(n,n,ksta:kend) :: vx,vy,vz,tmp1,tmp2
    REAL(KIND=GP)                                        :: dtfact
    REAL(KIND=GP),ALLOCATABLE  ,DIMENSION            (:) :: lid,gid

    ! ... x:
    CALL GPart_EulerToLag(this,this%lvx_,this%nparts_,vx,doupdate,tmp1,tmp2)
    ! ... y:
    CALL GPart_EulerToLag(this,this%lvy_,this%nparts_,vy,.false.,tmp1,tmp2)
    ! ... z:
    CALL GPart_EulerToLag(this,this%lvz_,this%nparts_,vz,.false.,tmp1,tmp2)

  END SUBROUTINE GPart_SetLagVec
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_EulerToLag(this,lag,nl,evar,doupdate,tmp1,tmp2)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : EulerToLag
!  DESCRIPTION: Computes the Eulerian to Lagrangian
!               transformation by interpolating Eulerian field
!               evar to current position of Lagrangian paricles in 
!               d.b. Array lag must be large enouch to accommodate 
!               max. no. particles; no checking is done. Note
!               that 'evar' array must have local dimensions 
!               for a real array in GHOST (n X n X (kend-ksta+1)).
!               Global computation of spline or other interpolation 
!               operator will be done here.
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    lag     : real array with the containing Lagrangian value
!              of evar field on output (OUT)
!    nl      : no. Lag. points in lag
!    evar    : Eulerian variable (IN)
!    doupdate: if true, do interp point update in interpolator; else don't
!    tmp1/2  : temp. arrays for interpolation
!    
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)                           :: this
    INTEGER      ,INTENT   (IN)                           :: nl
    LOGICAL      ,INTENT   (IN)                           :: doupdate
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(n,n,ksta:kend)  :: evar,tmp1,tmp2
    REAL(KIND=GP),INTENT(INOUT),DIMENSION           (nl)  :: lag
    INTEGER                                               :: j

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

  END SUBROUTINE GPart_EulerToLag
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_MakePeriodicP(this,px,py,pz,npdb,idir)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : MakePeriodicP
!  DESCRIPTION: Enforces periodic b.c.'s on particles in pdb
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    px      : particle x loc. d.b.
!    py      : particle x loc. d.b.
!    pz      : particle x loc. d.b.
!    npdb    : no. particles in pdb
!    idir    : first three bits provide which directions to periodize.
!              So, 1==>x, 2==>y, 4==>z; 3==>x & y; 7==>x & y & z
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)                      :: this
    INTEGER,INTENT(IN)                               :: idir,npdb
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)         :: px,py,pz

    INTEGER                                          :: j


    IF ( btest(idir,0) ) THEN
      DO j = 1, npdb
        px(j) = modulo(px(j)+2.0*this%gext_(1),this%gext_(1))
      ENDDO
    ENDIF
    
    IF ( btest(idir,1) ) THEN
      DO j = 1, npdb
        py(j) = modulo(py(j)+2.0*this%gext_(2),this%gext_(2))
      ENDDO
    ENDIF

    IF ( btest(idir,2) ) THEN
      DO j = 1, npdb
        pz(j) = modulo(pz(j)+2.0*this%gext_(3),this%gext_(3))
      ENDDO
    ENDIF

  END SUBROUTINE GPart_MakePeriodicP
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_MakePeriodicExt(this,v,nx,ny,kb,ke,nc)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : MakePeriodicExt
!  DESCRIPTION: Enforces periodic b.c.'s on extended field
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    v       : real field on extended grid 
!    nx,ny   : global size of field in x-y (including ghost zones)
!    kb,ke   : starting, ending z-indices of slab (including ghost zones)
!    nc      : index in x and y (and z) s.t. f(nc+1,:,:) = f(nx-nc,:,:), etc.
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)                      :: this
    INTEGER      ,INTENT(IN)                         :: nc,nx,ny,kb,ke
    REAL(KIND=GP),INTENT(INOUT)                      :: v(nx,ny,kb:ke)

    INTEGER                                          :: i,j,k

    ! Recall: nx, ny are the dimensions _including_ ghost zones:
    !
    ! Periodicity s.t.:
    !   | | [ | | | | ] | |
    !   a b       a b 
    DO k = kb,ke 
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

  END SUBROUTINE GPart_MakePeriodicExt
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_Delete(this,id,px,py,pz,npdb,nnew)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPart_Delete
!  DESCRIPTION: Removes from PDB NULL particles, concatenates list,
!               and sets new number of particles
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    id      : part ids
!    px,py pz: part. d.b.
!    npdb    : no. parts. in pdb
!    nnew    : no. non-NULL particles (set in GPartComm class)
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)                :: this
    INTEGER      ,INTENT   (IN)                :: npdb
    INTEGER      ,INTENT(INOUT),DIMENSION(npdb):: id
    INTEGER      ,INTENT  (OUT)                :: nnew
    INTEGER                                    :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(npdb):: px,py,pz

    nnew = 0
    DO i = 1, npdb
       IF ( this%id_(i) .NE. GPNULL ) nnew = nnew + 1
    ENDDO

    j     = 1
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

  END SUBROUTINE GPart_Delete
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_GetLocalWrk(this,id,lx,ly,lz,nl,gvdb,ngvdb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPart_GetLocalWrk
!  DESCRIPTION: Removes from PDB NULL particles, concatenates list,
!               and sets new number of particles
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    id      : part ids
!    lx,ly,lz: local part. d.b. vectors
!    nl      : no. parts. in local pdb
!    gvdb    : global VDB containing part. position records. Location
!              gives particle id.
!    ngvdb   : no. records in global VDB
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)                           :: this
    INTEGER      ,INTENT(INOUT)                           :: nl
    INTEGER      ,INTENT(INOUT),DIMENSION(this%maxparts_) :: id
    INTEGER      ,INTENT   (IN)                           :: ngvdb
    INTEGER                                               :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(this%maxparts_) :: lx,ly,lz
    REAL(KIND=GP),INTENT   (IN),DIMENSION(3,ngvdb)        :: gvdb

    nl = 0
    id = GPNULL
    DO j = 1, ngvdb
      IF ( gvdb(3,j).GE.this%lxbnds_(3,1) .AND. gvdb(3,j).LT.this%lxbnds_(3,2) ) THEN 
        nl = nl + 1
        id (nl) = j-1
        lx (nl) = gvdb(1,j)
        ly (nl) = gvdb(2,j)
        lz (nl) = gvdb(3,j)
      ENDIF
    ENDDO

  END SUBROUTINE GPart_GetLocalWrk
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPart_GetLocalWrk_aux(this,id,lx,ly,lz,tx,ty,tz,nl,gvdb,gtmp,ngvdb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPart_GetLocalWrk_aux
!  DESCRIPTION: Removes from PDB NULL particles, concatenates list,
!               and sets new number of particles. This auxiliary 
!               subroutines also updates arrays used during the 
!               intermediate steps of the RK solver, and is needed 
!               if local work is recomputed in the midde of a RK 
!               iteration.
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    id      : part ids
!    lx,ly,lz: local part. d.b. vectors
!    tx,ty,tz: local initial part. d.b. vectors
!    nl      : no. parts. in local pdb
!    gvdb    : global VDB containing part. position records. Location
!              gives particle id.
!    gtmp    : global VDB containing part. position records at the
!              beginning of the RK loop. Location gives particle id.
!    ngvdb   : no. records in global VDB
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)                           :: this
    INTEGER      ,INTENT(INOUT)                           :: nl
    INTEGER      ,INTENT(INOUT),DIMENSION(this%maxparts_) :: id
    INTEGER      ,INTENT   (IN)                           :: ngvdb
    INTEGER                                               :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(this%maxparts_) :: lx,ly,lz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(this%maxparts_) :: tx,ty,tz
    REAL(KIND=GP),INTENT   (IN),DIMENSION(3,ngvdb)        :: gvdb
    REAL(KIND=GP),INTENT   (IN),DIMENSION(3,ngvdb)        :: gtmp

    nl = 0
    id = GPNULL
    DO j = 1, ngvdb
      IF ( gvdb(3,j).GE.this%lxbnds_(3,1) .AND. gvdb(3,j).LT.this%lxbnds_(3,2) ) THEN 
        nl = nl + 1
        id (nl) = j-1
        lx (nl) = gvdb(1,j)
        ly (nl) = gvdb(2,j)
        lz (nl) = gvdb(3,j)
        tx (nl) = gtmp(1,j)
        ty (nl) = gtmp(2,j)
        tz (nl) = gtmp(3,j)
      ENDIF
    ENDDO

  END SUBROUTINE GPart_GetLocalWrk_aux
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPart_GetVDB(this,pdb,npdb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPart_GetVDB
!  DESCRIPTION: Gets particle d.b.
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    pdb     : part pdb, of size (3,npdb)
!    npdb    : size of pdb array (2nd dimension); must be >= maxparts_
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes

    IMPLICIT NONE 
    CLASS(GPart) ,INTENT(INOUT)                   :: this 
    INTEGER      ,INTENT   (IN)                   :: npdb
    INTEGER                                       :: j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(3,npdb):: pdb 

    IF ( this%iexchtype_.EQ.GPEXCHTYPE_NN ) THEN
      IF ( .NOT.GPart_PartNumConsistent(this,this%nparts_) ) THEN
          IF ( this%myrank_.eq.0 ) THEN
            WRITE(*,*) 'GPart_io_write_vel: Inconsistent particle count'
            STOP
        ENDIF
      ENDIF

      CALL this%gpcomm_%VDBSynch(pdb,this%maxparts_,this%id_, &
           this%px_,this%py_,this%pz_,this%nparts_,this%ptmp0_)
    ELSE
      DO j = 1, npdb
        pdb(1:3,j) = this%vdb_(1:3,j)
     ENDDO
    ENDIF

   END SUBROUTINE GPart_GetVDB
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPart_GetVel(this,lvel,nparts)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPart_GetVel
!  DESCRIPTION: Gets current particle velocities by doing 'synch'
!               of local velocities
!         
!  ARGUMENTS  :
!    this     : 'this' class instance (IN)
!    lvel     : part velocity array, of size (3,nparts)
!    nparts   : size of lvel array, must be >= maxparts_
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes

    IMPLICIT NONE 
    CLASS(GPart) ,INTENT(INOUT)                    :: this 
    INTEGER      ,INTENT   (IN)                    :: nparts
    INTEGER                                        :: j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(3,nparts):: lvel

    IF ( .NOT.GPart_PartNumConsistent(this,this%nparts_) ) THEN
      IF ( this%myrank_.eq.0 ) THEN
        WRITE(*,*) 'GPart_GetVel: Inconsistent particle count'
        STOP
      ENDIF
    ENDIF
    CALL this%gpcomm_%VDBSynch(lvel,this%maxparts_,this%id_, &
         this%lvx_,this%lvy_,this%lvz_,this%nparts_,this%ptmp0_)

   END SUBROUTINE GPart_GetVel
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  DOUBLE PRECISION FUNCTION GPart_GetTime(this,itime)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPart_GetTime
!  DESCRIPTION: gets elapsed time from timer index itime
!         
!  ARGUMENTS  :
!    this     : 'this' class instance (IN)
!    itime    : 'GPTIME' parameter (above)
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes

    IMPLICIT NONE
    CLASS(GPart) ,INTENT(INOUT)                    :: this
    INTEGER      ,INTENT   (IN)                    :: itime
    INTEGER                                        :: j

    IF ( itime.LT.GPTIME_STEP .OR. itime.GT.GPMAXTIMERS ) THEN
      WRITE(*,*)'GPart_GetTime: invalid time specification'
      STOP
    ENDIF

    GPart_GetTime = GTGetTime(this%htimers_(itime))

   END FUNCTION GPart_GetTime
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  INTEGER FUNCTION GPart_GetNParts(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPart_GetNParts
!  DESCRIPTION: Gets no. particles on grid
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!-----------------------------------------------------------------
     CLASS(GPart) ,INTENT(INOUT)                   :: this 
     INTEGER                                       :: ngp
    
     CALL MPI_ALLREDUCE(this%nparts_,ngp,1,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
     GPart_GetNParts = ngp

  END FUNCTION GPart_GetNParts
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  REAL FUNCTION GPart_GetLoadBal(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPart_GetLoadBal
!  DESCRIPTION: Gets current load (im)balance measure
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!-----------------------------------------------------------------
     CLASS(GPart) ,INTENT(INOUT)                   :: this 
     REAL                                          :: rbal      
     INTEGER                                       :: gnmax,gnmin
    
     CALL MPI_ALLREDUCE(this%nparts_,gnmin,1,MPI_INTEGER,MPI_MIN,this%comm_,this%ierr_)
     CALL MPI_ALLREDUCE(this%nparts_,gnmax,1,MPI_INTEGER,MPI_MAX,this%comm_,this%ierr_)
     rbal = real(gnmax) / (real(gnmin)+tiny(1.0))
     GPart_GetLoadBal = rbal

  END FUNCTION GPart_GetLoadBal
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  LOGICAL FUNCTION GPart_PartNumConsistent(this,nlocal)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPart_PartNumConsistent
!  DESCRIPTION: Checks that sum of local particle counts equals
!               maxparts_
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!-----------------------------------------------------------------
     CLASS(GPart) ,INTENT(INOUT)                   :: this 
     REAL                                          :: rbal      
     INTEGER,INTENT(IN)                            :: nlocal
     INTEGER                                       :: ng
    
     CALL MPI_ALLREDUCE(nlocal,ng,1,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
     GPart_PartNumConsistent = ng .EQ. this%maxparts_

  END FUNCTION GPart_PartNumConsistent
!-----------------------------------------------------------------
!-----------------------------------------------------------------

! SUBROUTINE GPart_GetParts(this,)
!!-----------------------------------------------------------------
!!-----------------------------------------------------------------
!!  METHOD     : GPart_GetParts
!!  DESCRIPTION: Gets particle d.b.
!!  ARGUMENTS  :
!!    this    : 'this' class instance (IN)
!!    id      : part ids
!!    lx,ly,lz: local part. d.b. vectors
!!    nl      : no. parts. in local pdb
!!    gvdb    : global VDB containing part. position records. Location
!!              gives particle id.
!!    ngvdb   : no. records in global VDB
!!-----------------------------------------------------------------
!    USE fprecision
!    USE commtypes
!
!    IMPLICIT NONE 
!    CLASS(GPart) ,INTENT(INOUT)                   :: this 
!    INTEGER      ,INTENT  (OUT)                   :: nl
!    INTEGER      ,INTENT  (OUT),DIMENSION(nl)     :: id
!    INTEGER      ,INTENT   (IN)                   :: ngvdb
!    INTEGER                                       :: i,j
!    REAL(KIND=GP),INTENT(OUT),DIMENSION(nl)       :: lx,ly,lz
!    REAL(KIND=GP),INTENT (IN),DIMENSION(3,ngvdb)  :: gvdb 
!
!
!  END SUBROUTINE GPart_GetParts
!-----------------------------------------------------------------
!-----------------------------------------------------------------


END MODULE class_GPart
