!=================================================================
! GHOST GPartComm particles communication class. It handles
!       two types of exchanges: the particles, and the
!       velocity data used to update the particle positions,
!       (the 'ghost' zone data) and separate interfaces are
!       provided for each. The velocity data can be exchanged
!       for _any_ other field variable as well, although the
!       'multi-field' interfaces may not be appropriate for it.
!
! 2013 D. Rosenberg
!      ORNL: NCCS
!
! 15 Aug 2011: Initial version
!=================================================================
MODULE class_GPartComm
      USE fprecision
      USE commtypes
      USE gtimer
      IMPLICIT NONE

      INTEGER,PARAMETER,PUBLIC                       :: GPNULL=-1          ! particle NULL value
      INTEGER,PARAMETER,PUBLIC                       :: GPCOMM_INTRFC_SF=0 ! single-field interface
      INTEGER,PARAMETER,PUBLIC                       :: GPCOMM_INTRFC_MF=1 ! multi-field interface
      INTEGER,PARAMETER,PUBLIC                       :: GPEXCH_INIT = 0    ! starts particle exchange
      INTEGER,PARAMETER,PUBLIC                       :: GPEXCH_UPDT = 1    ! continues part. exchange
      INTEGER,PARAMETER,PUBLIC                       :: GPEXCH_END  = 2    ! finishes particle exchange
      INTEGER,PARAMETER,PUBLIC                       :: GPEXCH_UNIQ = 3    ! exchange only positions
!$    INTEGER,PARAMETER,PUBLIC                       :: NMIN_OMP    = 10000! min. no. of op. for threadingi
      INTEGER,PUBLIC                                 :: MPI_GPDataType=-1,MPI_GPDataPackType=-1
!!!
      INTEGER                                        :: UNPACK_REP = 0, UNPACK_SUM = 1
!!!
      PRIVATE
      TYPE, BIND(C), PUBLIC :: GPData
        INTEGER        :: idp_
        REAL(KIND=GP)  :: rp_(3)
      END TYPE GPData

      TYPE, PUBLIC :: GPartComm
        PRIVATE
        ! Member data:
        INTEGER                                      :: intrfc_   ! if >=1 use multi-field interface
        INTEGER                                      :: maxparts_  ,nbuff_     ,nd_(3)   ,nzghost_
        INTEGER                                      :: nbsnd_     ,ntsnd_     ,nbrcv_   ,ntrcv_
        INTEGER                                      :: nprocs_    ,myrank_    ,comm_
        INTEGER                                      :: csize_     ,nstrip_
        INTEGER                                      :: ntop_      ,nbot_      ,ierr_
        INTEGER                                      :: iextperp_  ,ksta_      ,kend_    ,nth_
        INTEGER                                      :: hcomm_
        LOGICAL                                      :: btransinit_
        INTEGER, ALLOCATABLE, DIMENSION(:,:)         :: ibsnd_     ,itsnd_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: ibrcv_     ,itrcv_     ,nbbrcv_  ,ntbrcv_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: ibsh_      ,itsh_      ,ibrh_    ,itrh_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: igsh_      ,igrh_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: itypekp_   ,itypeip_   ,itypea_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: ibsndnz_   ,itsndnz_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: itop_      ,ibot_
        INTEGER, ALLOCATABLE, DIMENSION(:,:)         :: ibsnddst_  ,itsnddst_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: ibrcvnz_   ,itrcvnz_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: ibsndp_    ,ibrcvp_    ,itsndp_  ,itrcvp_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: oldid_
        INTEGER, DIMENSION (MPI_STATUS_SIZE)         :: istatus_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:)   :: sbbuff_    ,stbuff_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:)   :: rbbuff_    ,rtbuff_
        TYPE(GPData) , ALLOCATABLE, DIMENSION(:)     :: sbbuffp_   ,stbuffp_
        TYPE(GPData) , ALLOCATABLE, DIMENSION(:)     :: rbbuffp_   ,rtbuffp_
!!!
        INTEGER                                      :: nbret_     ,ntret_
        INTEGER, ALLOCATABLE, DIMENSION(:)           :: ibretp_    ,itretp_    ,ibretnz_ ,itretnz_
        INTEGER, ALLOCATABLE, DIMENSION(:,:)         :: ibret_     ,itret_     ,ibretdst_,itretdst_
!!!
      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: GPartComm_ctor
        FINAL            :: GPartComm_dtor
        PROCEDURE,PUBLIC :: Init              => GPartComm_Init
        PROCEDURE,PUBLIC :: GetNumGhost       => GPartComm_GetNumGhost
        PROCEDURE,PUBLIC :: GTranspose        => GPartComm_Transpose
        PROCEDURE,PUBLIC :: GITranspose       => GPartComm_ITranspose
        PROCEDURE,PUBLIC :: VDBSynch          => GPartComm_VDBSynch
        PROCEDURE,PUBLIC :: VDBSynch_t0       => GPartComm_VDBSynch_t0
        PROCEDURE,PUBLIC :: LagSynch          => GPartComm_LagSynch
        PROCEDURE,PUBLIC :: LagSynch_t0       => GPartComm_LagSynch_t0
        PROCEDURE,PUBLIC :: SetCacheParam     => GPartComm_SetCacheParam
        PROCEDURE,PUBLIC :: PartExchangePDB   => GPartComm_PartExchangePDB
        PROCEDURE,PUBLIC :: IdentifyTaskV     => GPartComm_IdentifyTaskV
        PROCEDURE,PUBLIC :: PartScatterV      => GPartComm_PartScatterV
        PROCEDURE,PUBLIC :: PartExchangeV     => GPartComm_PartExchangeV
        PROCEDURE,PUBLIC :: IdentifyExchV     => GPartComm_IdentifyExchV
        PROCEDURE,PUBLIC :: SlabDataExchangeMF=> GPartComm_SlabDataExchangeMF
        PROCEDURE,PUBLIC :: SlabDataExchangeSF=> GPartComm_SlabDataExchangeSF
        PROCEDURE,PUBLIC :: ConcatPDB         => GPartComm_ConcatPDB
        PROCEDURE,PUBLIC :: ConcatV           => GPartComm_ConcatV
        PROCEDURE,PUBLIC :: Copy2Ext          => GPartComm_Copy2Ext
        GENERIC  ,PUBLIC :: Concat            => ConcatPDB,ConcatV
        PROCEDURE,PUBLIC :: ResizeArrays      => GPartComm_ResizeArrays
!!!
        PROCEDURE,PUBLIC :: SlabDataReturnMF  => GPartComm_SlabDataReturnMF
        PROCEDURE,PUBLIC :: SlabDataReturnSF  => GPartComm_SlabDataReturnSF
        PROCEDURE,PUBLIC :: Copy2Reg          => GPartComm_Copy2Reg
        PROCEDURE,PUBLIC :: AllocRetArrays    => GPartComm_AllocRet
!!!
      END TYPE GPartComm

      PUBLIC :: Resize_ArrayRank1,Resize_ArrayRank2
      PUBLIC :: Resize_IntArray,Resize_IntArrayRank2
      PUBLIC :: Resize_ArrayRank2Transposed

!      INCLUDE 'gfieldcomm_dtype.f90'

      PRIVATE :: GPartComm_Init              , GPartComm_AllocRet
      PRIVATE :: GPartComm_SlabDataExchangeMF, GPartComm_SlabDataExchangeSF
      PRIVATE :: GPartComm_LocalDataExchMF   , GPartComm_LocalDataExchSF
      PRIVATE :: GPartComm_PartExchangePDB   , GPartComm_PartExchangeV
      PRIVATE :: GPartComm_Transpose         , GPartComm_GetNumGhost
      PRIVATE :: GPartComm_PackMF            , GPartComm_UnpackMF
      PRIVATE :: GPartComm_PackSF            , GPartComm_UnpackSF
      PRIVATE :: GPartComm_PPackPDB          , GPartComm_PUnpackPDB
      PRIVATE :: GPartComm_PPackV            , GPartComm_PUnpackV
      PRIVATE :: GPartComm_SetCacheParam
!!!
      PRIVATE :: GPartComm_PackRetSF         , GPartComm_PackRetMF
      PRIVATE :: GPartComm_UnpackRetSF       , GPartComm_UnpackRetMF
      PRIVATE :: GPartComm_LocalDataRetMF    , GPartComm_LocalDataRetSF
!!!      
!      INCLUDE 'gfieldcomm_private.f90'

! Methods:
  CONTAINS

  SUBROUTINE GPartComm_ctor(this,intrface,maxparts,nd,nzghost,comm,hcomm)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Main explicit constructor
!  ARGUMENTS:
!    this    : 'this' class instance
!    intrface: which interface (MF (>=1) or SF (0) is expected. MF uses more
!              memory, but should be faster.
!    nparts  : no. particles allowed on grid.
!    nd(3)   : x- ,y- , and z- (global) dimensions of data
!    nzghost : 'z' : no. slices of each slab required to
!              build 'ghost' zones.  If there are fewer slices on
!              adjacent tasks, method will go 'adjacent' tasks to find
!              the information to fill ghost zones.
!    comm    : MPI communicator
!    hcomm   : externally-managed comm-timer handle; must be non-null on entry
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT):: this
    INTEGER, INTENT(IN)           :: intrface,maxparts,nd(3),nzghost
    INTEGER, INTENT(IN)           :: comm,hcomm
!$  INTEGER, EXTERNAL             :: omp_get_max_threads

    this%intrfc_    = intrface
    this%maxparts_  = maxparts
    this%nd_        = nd
    this%nzghost_   = nzghost
    this%comm_      = comm;
    this%csize_     = 8;
    this%nstrip_    = 1;
    this%iextperp_  = 0;     ! set extended grid in perp direction (x-y) too?
    this%nth_       = 1
    IF ( GTValidHandle(hcomm).NE.GTERR_GOOD_HANDLE ) THEN
      WRITE(*,*) 'GPPartComm_ctor: invalid comm timer handle: ',hcomm
      STOP
    ENDIF
    this%hcomm_     = hcomm;
!$  this%nth_       = omp_get_max_threads()

    CALL MPI_COMM_SIZE(this%comm_,this%nprocs_,this%ierr_)
    IF ( this%ierr_ .NE. MPI_SUCCESS ) THEN
      WRITE(*,*) 'GPartComm::ctor: MPI_COMM_SIZE: err=',this%ierr_
      STOP
    ENDIF
    CALL MPI_COMM_RANK(this%comm_,this%myrank_,this%ierr_)
    IF ( this%ierr_ .NE. MPI_SUCCESS ) THEN
      WRITE(*,*) 'GPartComm::ctor: MPI_COMM_RANK: err=',this%ierr_
      STOP
    ENDIF
    this%btransinit_ = .FALSE.

    IF ( this%intrfc_ .GE. 1 ) THEN
!      this%nbuff_  = MAX(maxparts*4,3*nd(1)*nd(2)*nzghost+nzghost+1) + 1
      this%nbuff_  = 3*nd(1)*nd(2)*nzghost+nzghost+2
    ELSE
!      this%nbuff_  = MAX(maxparts*4,  nd(1)*nd(2)*nzghost+nzghost+1) + 1
      this%nbuff_  =   nd(1)*nd(2)*nzghost+nzghost+2
    ENDIF

    CALL GPartComm_Init(this)

  END SUBROUTINE GPartComm_ctor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

 SUBROUTINE GPartComm_AllocRet(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Allocate return arrays for field communicator. Should be called
!  after calling GPartComm_Init.
!
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPartComm), INTENT(INOUT)       :: this
    INTEGER                               :: jf,kg,k,n2p,nt
    INTEGER,ALLOCATABLE,DIMENSION(:)      :: jfwd,kfend,kfsta,nzf

    n2p = this%nd_(3)/this%nprocs_
    nt  = (this%nzghost_+n2p-1)/n2p  ! max no. tasks needed for ghost zones

    ALLOCATE(jfwd (nt))
    ALLOCATE(kfend(nt))
    ALLOCATE(kfsta(nt))
    ALLOCATE(nzf  (nt))

    ALLOCATE(this%ibretp_(nt))
    ALLOCATE(this%itretp_(nt))
    ALLOCATE(this%ibretnz_(nt))
    ALLOCATE(this%itretnz_(nt))
    ALLOCATE(this%ibret_(nt,this%nzghost_+1))
    ALLOCATE(this%itret_(nt,this%nzghost_+1))
    ALLOCATE(this%ibretdst_(nt,this%nzghost_+1))
    ALLOCATE(this%itretdst_(nt,this%nzghost_+1))

    ! *** Find top neighbors to return to:
    DO jf = 1, nt
      jfwd(jf) = modulo(this%myrank_+jf,this%nprocs_)
      CALL range(1,this%nd_(3),this%nprocs_,jfwd(jf),kfsta(jf),kfend(jf))
      nzf (jf) = kfend(jf) - kfsta(jf) + 1
    ENDDO

    this%ntret_ = 1
    this%itretp_(this%ntret_) = jfwd(this%ntret_)
    this%itretnz_(this%ntret_) = 0
    k = 1
    DO kg = 1, this%nzghost_
      IF ( k .GT. nzf(this%ntret_) ) THEN
        k = k - nzf(this%ntret_)
        this%ntret_ = this%ntret_ + 1
        this%itretp_(this%ntret_)  = jfwd(this%ntret_)
        this%itretnz_(this%ntret_) = 0
      ENDIF
      ! local z-index to be returned to top (in _extended_ grid)
      this%itret_   (this%ntret_,k) = this%nzghost_+this%kend_-this%ksta_+1+kg
      ! Destination z-indices in _regular_ grid for top return.
      ! These indices should be in local--not global--form:
      this%itretdst_(this%ntret_,k) = k
      ! Set no. z-indices to return to top task
      this%itretnz_ (this%ntret_  ) = this%itretnz_(this%ntret_) + 1
      k = k + 1
    ENDDO

    ! *** Find bottom neighbors to return to:
    DO jf = 1, nt
      jfwd(jf) = modulo(this%myrank_-jf+this%nprocs_,this%nprocs_)
      CALL range(1,this%nd_(3),this%nprocs_,jfwd(jf),kfsta(jf),kfend(jf))
      nzf (jf) = kfend(jf) - kfsta(jf) + 1
    ENDDO

    this%nbret_ = 1
    this%ibretp_(this%nbret_) = jfwd(this%nbret_)
    this%ibretnz_(this%nbret_) = 0
    k = 1
    DO kg = 1, this%nzghost_
      IF ( k .GT. nzf(this%nbret_) ) THEN
        k = k - nzf(this%nbret_)
        this%nbret_ = this%nbret_ + 1
        this%ibretp_(this%nbret_)  = jfwd(this%nbret_)
        this%ibretnz_(this%nbret_) = 0
      ENDIF
      ! local z-index to be returned to top (in _extended_ grid)
      this%ibret_   (this%nbret_,k) = this%nzghost_ - kg + 1
      ! Destination z-indices in _regular_ grid for bottom return.
      ! These indices should be in local--not global--form:
      this%ibretdst_(this%nbret_,k) = nzf(this%nbret_) + 1 - k
      ! Set no. z-indices to return to bottom task
      this%ibretnz_ (this%nbret_  ) = this%ibretnz_(this%nbret_) + 1
      k = k + 1
    ENDDO

    DEALLOCATE(jfwd,kfend,kfsta,nzf)

  END SUBROUTINE GPartComm_AllocRet
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_dtor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : dtor
!  DESCRIPTION: Destructor
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    TYPE(GPartComm),INTENT(INOUT)        :: this

    CALL GPartComm_DoDealloc(this)

  END SUBROUTINE GPartComm_dtor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_DoDealloc(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : DoDealloc
!  DESCRIPTION: Does allocation of class member data
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    TYPE(GPartComm),INTENT(INOUT)        :: this

    IF ( ALLOCATED   (this%sbbuff_) ) DEALLOCATE  (this%sbbuff_)
    IF ( ALLOCATED   (this%stbuff_) ) DEALLOCATE  (this%stbuff_)
    IF ( ALLOCATED   (this%rbbuff_) ) DEALLOCATE  (this%rbbuff_)
    IF ( ALLOCATED   (this%rtbuff_) ) DEALLOCATE  (this%rtbuff_)
    IF ( ALLOCATED  (this%sbbuffp_) ) DEALLOCATE (this%sbbuffp_)
    IF ( ALLOCATED  (this%stbuffp_) ) DEALLOCATE (this%stbuffp_)
    IF ( ALLOCATED  (this%rbbuffp_) ) DEALLOCATE (this%rbbuffp_)
    IF ( ALLOCATED  (this%rtbuffp_) ) DEALLOCATE (this%rtbuffp_)
    IF ( ALLOCATED     (this%itop_) ) DEALLOCATE    (this%itop_)
    IF ( ALLOCATED     (this%ibot_) ) DEALLOCATE    (this%ibot_)

    IF ( ALLOCATED   (this%ibrcvp_) ) DEALLOCATE  (this%ibrcvp_)
    IF ( ALLOCATED   (this%ibsndp_) ) DEALLOCATE  (this%ibsndp_)
    IF ( ALLOCATED   (this%itrcvp_) ) DEALLOCATE  (this%itrcvp_)
    IF ( ALLOCATED   (this%itsndp_) ) DEALLOCATE  (this%itsndp_)
    IF ( ALLOCATED    (this%ibrcv_) ) DEALLOCATE   (this%ibrcv_)
    IF ( ALLOCATED   (this%nbbrcv_) ) DEALLOCATE  (this%nbbrcv_)
    IF ( ALLOCATED    (this%itrcv_) ) DEALLOCATE   (this%itrcv_)
    IF ( ALLOCATED   (this%ntbrcv_) ) DEALLOCATE  (this%ntbrcv_)
    IF ( ALLOCATED    (this%ibsnd_) ) DEALLOCATE   (this%ibsnd_)
    IF ( ALLOCATED    (this%itsnd_) ) DEALLOCATE   (this%itsnd_)
    IF ( ALLOCATED     (this%ibrh_) ) DEALLOCATE    (this%ibrh_)
    IF ( ALLOCATED     (this%itrh_) ) DEALLOCATE    (this%itrh_)
    IF ( ALLOCATED     (this%ibsh_) ) DEALLOCATE    (this%ibsh_)
    IF ( ALLOCATED     (this%itsh_) ) DEALLOCATE    (this%itsh_)
    IF ( ALLOCATED     (this%igrh_) ) DEALLOCATE    (this%igrh_)
    IF ( ALLOCATED     (this%igsh_) ) DEALLOCATE    (this%igsh_)
    IF ( ALLOCATED  (this%itypekp_) ) DEALLOCATE (this%itypekp_)
    IF ( ALLOCATED  (this%itypeip_) ) DEALLOCATE (this%itypeip_)
    IF ( ALLOCATED   (this%itypea_) ) DEALLOCATE  (this%itypea_)
    IF ( ALLOCATED (this%ibsnddst_) ) DEALLOCATE(this%ibsnddst_)
    IF ( ALLOCATED (this%itsnddst_) ) DEALLOCATE(this%itsnddst_)
    IF ( ALLOCATED  (this%ibrcvnz_) ) DEALLOCATE (this%ibrcvnz_)
    IF ( ALLOCATED  (this%itrcvnz_) ) DEALLOCATE (this%itrcvnz_)
    IF ( ALLOCATED  (this%ibsndnz_) ) DEALLOCATE (this%ibsndnz_)
    IF ( ALLOCATED  (this%itsndnz_) ) DEALLOCATE (this%itsndnz_)

    IF ( ALLOCATED  (this%oldid_)   ) DEALLOCATE   (this%oldid_)
!!!
    IF ( ALLOCATED (this%ibretp_  ) ) DEALLOCATE (this%ibretp_  )
    IF ( ALLOCATED (this%itretp_  ) ) DEALLOCATE (this%itretp_  )
    IF ( ALLOCATED (this%ibretnz_ ) ) DEALLOCATE (this%ibretnz_ )
    IF ( ALLOCATED (this%itretnz_ ) ) DEALLOCATE (this%itretnz_ )
    IF ( ALLOCATED (this%ibret_   ) ) DEALLOCATE (this%ibret_   )
    IF ( ALLOCATED (this%itret_   ) ) DEALLOCATE (this%itret_   )
    IF ( ALLOCATED (this%ibretdst_) ) DEALLOCATE (this%ibretdst_)
    IF ( ALLOCATED (this%itretdst_) ) DEALLOCATE (this%itretdst_)
!!!

  END SUBROUTINE GPartComm_DoDealloc
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_Init(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Init
!  DESCRIPTION: Initializes particle locations before integration.
!               Call after construction.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    comm    : MP:I communicator
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)        :: this
    INTEGER                               :: i,ibrank,itrank,j,jf,k,kf,kend,ksta
    INTEGER                               :: kbend,kbsta,ktend,ktsta
    INTEGER                               :: jbr,jtr,n2p,nt,nxy
    INTEGER,ALLOCATABLE,DIMENSION(:)      :: jfwd,kfend,kfsta,nzf

    ! Compute the no. sendto and recv from tasks there are:

    ! If there aren't enough 'slices' with nearest neighbors to
    ! fill ghost zones, go to next furthest tasks, etc., to fill:
    this%nbsnd_ = 0
    this%ntsnd_ = 0
    this%nbrcv_ = 0
    this%ntrcv_ = 0
    n2p = this%nd_(3)/this%nprocs_
    nt  = (this%nzghost_+n2p-1)/n2p  ! max no. tasks needed for ghost zones

    CALL GPartComm_DoDealloc(this)

    ALLOCATE(this%sbbuff_ (this%nbuff_,nt))
    ALLOCATE(this%stbuff_ (this%nbuff_,nt))
    ALLOCATE(this%rbbuff_ (this%nbuff_,nt))
    ALLOCATE(this%rtbuff_ (this%nbuff_,nt))
    ALLOCATE(this%sbbuffp_(this%maxparts_))
    ALLOCATE(this%stbuffp_(this%maxparts_))
    ALLOCATE(this%rbbuffp_(this%maxparts_))
    ALLOCATE(this%rtbuffp_(this%maxparts_))

    ALLOCATE(this%ibot_(this%maxparts_))
    ALLOCATE(this%itop_(this%maxparts_))
    ALLOCATE(this%oldid_(this%maxparts_))

    ALLOCATE(this%ibrcvp_(nt))
    ALLOCATE(this%ibsndp_(nt))
    ALLOCATE(this%itrcvp_(nt))
    ALLOCATE(this%itsndp_(nt))
    ALLOCATE(this%ibrcv_(nt))
    ALLOCATE(this%nbbrcv_(nt))
    ALLOCATE(this%itrcv_(nt))
    ALLOCATE(this%ntbrcv_(nt))
    ALLOCATE(this%ibsnd_(nt,this%nzghost_+1))
    ALLOCATE(this%itsnd_(nt,this%nzghost_+1))
    ALLOCATE(this%ibrh_(nt))
    ALLOCATE(this%itrh_(nt))
    ALLOCATE(this%ibsh_(nt))
    ALLOCATE(this%itsh_(nt))
    ALLOCATE(this%ibrcvnz_(nt))
    ALLOCATE(this%itrcvnz_(nt))
    ALLOCATE(this%ibsndnz_(nt))
    ALLOCATE(this%itsndnz_(nt))
    ALLOCATE(this%igrh_(0:this%nprocs_-1))
    ALLOCATE(this%igsh_(0:this%nprocs_-1))
    ALLOCATE(this%itypekp_(0:this%nprocs_-1))
    ALLOCATE(this%itypeip_(0:this%nprocs_-1))
    ALLOCATE(this%itypea_(0:this%nprocs_-1))
    ALLOCATE(this%ibsnddst_(nt,this%nzghost_+1))
    ALLOCATE(this%itsnddst_(nt,this%nzghost_+1))

    ! Initialize all task/neighbor  lists with GPNULL:
    this%ibrcv_   =GPNULL; this%itrcv_   =GPNULL; this%ibsnd_  =GPNULL; this%itsnd_ =GPNULL;
    this%ibsnddst_=GPNULL; this%itsnddst_=GPNULL;

    ! Get global z-bounds on this rank:
    CALL range(1,this%nd_(3),this%nprocs_,this%myrank_,ksta,kend)
    this%ksta_ = ksta
    this%kend_ = kend

    this%ntsnd_   = 0
    this%ntrcv_   = 0
    this%nbrcv_   = 0
    this%nbsnd_   = 0

    this%itrcvnz_ = 0
    this%itsndnz_ = 0
    this%ibrcvnz_ = 0
    this%ibsndnz_ = 0

    ALLOCATE(jfwd (nt))
    ALLOCATE(kfend(nt))
    ALLOCATE(kfsta(nt))
    ALLOCATE(nzf  (nt))

    ! *** Find top neighbors to send to:
    DO jf = 1, nt
      jfwd(jf) = modulo(this%myrank_+jf,this%nprocs_)
      CALL range(1,this%nd_(3),this%nprocs_,jfwd(jf),kfsta(jf),kfend(jf))
      nzf (jf) = kfend(jf)-kfsta(jf)+1
    ENDDO

    ! Send to myrank + 1 at top:
    this%itsndp_(1) = jfwd(1)
    DO k = 1, min(kend-ksta+1, this%nzghost_)
      this%itsnd_   (1,k) = kend-ksta+2-k ! local z-index to be sent to top
      ! Destination z-indices in _extended_ grid for itsnd.
      ! These indices should be in local--not global--form:
      this%itsnddst_(1,k) = this%nzghost_-k+1
      ! Set no. z-indices to send to top task; used to
      ! set where in recv buffer to store recvd data:
      this%itsndnz_   (1) = this%itsndnz_(1) + 1
    ENDDO
    this%ntsnd_ = 1

    ! Send to other top neighbors that need data:
    DO jf = 2, nt
      IF ( nzf(jf-1) .LT. this%nzghost_ ) THEN
        kf = this%nzghost_-nzf(jf-1)
        IF ( kf .GE. 1 ) THEN
          this%ntsnd_ =  this%ntsnd_ + 1
          this%itsndp_(this%ntsnd_) = jfwd(jf)   ! top task to send to
          DO k = 1, min(kend-ksta+1, kf)
            this%itsnd_   (this%ntsnd_,k) = kend-ksta+2-k ! local z-index to be sent to top
            this%itsnddst_(this%ntsnd_,k) = this%nzghost_-kf-k+1 ! local destination z-index
            this%itsndnz_ (this%ntsnd_  ) = this%itsndnz_(this%ntsnd_) + 1 ! gives position in top recv buffer
          ENDDO
        ENDIF
      ENDIF
    ENDDO

    ! Find top neighbors to receive from:
    jf = 0
    DO WHILE ( this%ntrcv_.LE.nt .AND. jf.LT.this%nzghost_ )
      itrank = modulo(this%myrank_+this%ntrcv_+1,this%nprocs_)
      CALL range(1,this%nd_(3),this%nprocs_,itrank,ktsta,ktend)
      this%itrcvp_(this%ntrcv_+1) = itrank    !top rcv task
      k = 1
      DO WHILE ( k.LE.ktend-ktsta+1 .AND. jf.LT.this%nzghost_ )
        jf = jf + 1
        this%itrcvnz_  (this%ntrcv_+1)    = this%itrcvnz_(this%ntrcv_+1) + 1
        k = k + 1
      ENDDO
      this%ntrcv_ = this%ntrcv_ + 1 ! no. procs to recv from at top
    ENDDO
    IF ( jf.NE.this%nzghost_ ) THEN
       WRITE(*,*) 'GPartComm_Init: top neighbor data incompatible with interpolation order'
       WRITE(*,*) 'GPartComm_Init: nghost=',this%nzghost_, ' nt=',jf
       STOP
    ENDIF

    ! *** Find bottom neighbors to send to:
    DO jf = 1, nt
      jfwd(jf) = modulo(this%myrank_-jf+this%nprocs_,this%nprocs_)
      CALL range(1,this%nd_(3),this%nprocs_,jfwd(jf),kfsta(jf),kfend(jf))
      nzf (jf) = kfend(jf)-kfsta(jf)+1
    ENDDO

    ! Send to myrank - 1 at bottom :
    this%ibsndp_(1) = jfwd(1)
    DO k = 1, min(kend-ksta+1, this%nzghost_)
      this%ibsnd_   (1,k) = k ! local z-index to be sent to bottom
      ! Destination z-indices in _extended_ grid for bottom snd.
      ! These indices should be in local--not global--form:
      this%ibsnddst_(1,k) = this%nzghost_+nzf(1)+k
      ! Set no. z-indices to send to bottom task; used to
      ! set where in recv buffer to store recvd data:
      this%ibsndnz_   (1) = this%ibsndnz_(1) + 1
    ENDDO
    this%nbsnd_ = 1

    ! Send to other bottom neighbors that need data:
    DO jf = 2, nt
      IF ( nzf(jf-1) .LT. this%nzghost_ ) THEN
        kf = this%nzghost_-nzf(jf-1)
        IF ( kf .GE. 1 ) THEN
          this%nbsnd_ =  this%nbsnd_ + 1
          this%ibsndp_(this%nbsnd_) = jfwd(jf)   ! bottom task to send to
          DO k = 1, min(kend-ksta+1, kf)
            this%ibsnd_   (this%nbsnd_,k) = k ! local z-index to be sent to bottom
            this%ibsnddst_(this%nbsnd_,k) = this%nzghost_+nzf(jf)+kf+k ! local destination z-index
            this%ibsndnz_ (this%nbsnd_  ) = this%ibsndnz_(this%nbsnd_) + 1 ! gives position in bottom recv buffer
          ENDDO
        ENDIF
      ENDIF
    ENDDO

    ! Find bottom receives:
    jf = 0
    DO WHILE ( this%nbrcv_.LE.nt .AND. jf.LT.this%nzghost_ )
      ibrank = modulo(this%myrank_-this%nbrcv_-1+this%nprocs_,this%nprocs_)
      CALL range(1,this%nd_(3),this%nprocs_,ibrank,kbsta,kbend)
      this%ibrcvp_(this%nbrcv_+1) = ibrank    !bottom rcv task
      k = 1
      DO WHILE ( k.LE.kbend-kbsta+1 .AND. jf.LT.this%nzghost_ )
        jf = jf + 1
        this%ibrcvnz_  (this%nbrcv_+1)    = this%ibrcvnz_(this%nbrcv_+1) + 1
        k = k + 1
      ENDDO
      this%nbrcv_ = this%nbrcv_ + 1 ! no. procs to recv from at bottom
    ENDDO
    IF ( jf.NE.this%nzghost_ ) THEN
       WRITE(*,*) 'GPartComm_Init: bottom neighbor data incompatible with interpolation order'
       WRITE(*,*) 'GPartComm_Init: nghost=',this%nzghost_, ' nb=',jf
       STOP
    ENDIF

    ! Indices in recv buff to put data recvd from task j;
    ! includes 2 integer header:
!
    ! For multifield interfaces, the rcv buff starting indices are different:
    ! than for single field interface:
    nxy = this%nd_(1)*this%nd_(2)
    IF ( this%intrfc_ .GT. 0 ) THEN  ! multi-field interface
      DO j=1,nt
        jtr = 0; jbr = 0
        DO i = 1,j-1
          jtr = jtr + this%itrcvnz_(j)
          jbr = jbr + this%ibrcvnz_(j)
        ENDDO
        this%itrcv_  (j)    = jtr * (3*nxy+1)+1
        this%ibrcv_  (j)    = jbr * (3*nxy+1)+1
        this%ntbrcv_ (j)    = this%itrcvnz_(j) * (3*nxy+1)+1
        this%nbbrcv_ (j)    = this%ibrcvnz_(j) * (3*nxy+1)+1
      ENDDO
    ELSE                              ! single-field interface
      DO j=1,nt
        jtr = 0; jbr = 0
        DO i = 1,j-1
          jtr = jtr + this%itrcvnz_(i)
          jbr = jbr + this%ibrcvnz_(i)
        ENDDO
        ! where to put the data from other procs in this rcv buffer:
        this%itrcv_  (j)    = jtr * (nxy+1)+1
        this%ibrcv_  (j)    = jbr * (nxy+1)+1
        this%ntbrcv_ (j)    = this%itrcvnz_(j) * (nxy+1)+1
        this%nbbrcv_ (j)    = this%ibrcvnz_(j) * (nxy+1)+1
      ENDDO
    ENDIF

    DEALLOCATE(jfwd,kfend,kfsta,nzF)
   
    CALL GPartComm_InitMPITypes()

  END SUBROUTINE GPartComm_Init
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_InitMPITypes()
!-----------------------------------------------------------------
!  METHOD     : InitMPITypes
!  DESCRIPTION: Creates MPI type equivalent for GPData and its
!               corresponding array version.
!  ARGUMENTS  :
!    this              : 'this' class instance (IN)
!-----------------------------------------------------------------
    TYPE(GPData)                   :: gpdat,gpdatarr(2)
    INTEGER                        :: blocklen(2),types(2),ierr
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp(2)

    blocklen(1) = 1
    blocklen(2) = 3
    types(1) = MPI_INTEGER
    types(2) = GC_REAL

    CALL MPI_GET_ADDRESS(gpdat%idp_, disp(1), ierr)
    CALL MPI_GET_ADDRESS(gpdat%rp_ , disp(2), ierr)
    disp(2) = disp(2) - disp(1)
    disp(1) = 0
    CALL MPI_TYPE_CREATE_STRUCT(2,blocklen,disp,types,MPI_GPDataType,ierr) 
    CALL MPI_TYPE_COMMIT(MPI_GPDataType, ierr)

    CALL MPI_GET_ADDRESS(gpdatarr(1), disp(1), ierr) 
    CALL MPI_GET_ADDRESS(gpdatarr(2), disp(2), ierr) 
    CALL MPI_TYPE_CREATE_RESIZED(MPI_GPDataType,0,disp(2)-disp(1), &
                                 MPI_GPDataPackType,ierr) 
    CALL MPI_TYPE_COMMIT(MPI_GPDataPackType, ierr) 

  END SUBROUTINE
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_LocalDataExchMF(this,vxext,vyext,vzext,vx,vy,vz)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : LocalDataExch
!  DESCRIPTION: Does 'bdy exchange' of velocity component, when there's
!               only a single MPI task.
!  ARGUMENTS  :
!    this              : 'this' class instance (IN)
!    vxext,vyext,vzext : Eulerian velocity components returned on extended
!                        grid (that used to hold ghost zones). Only z-conditions
!                        are imposed; lateral periodicity is not handled here.
!                        Lateral ghost zones can be accounted for by setting
!                        this%iextperp_=1 in contructor.
!    vx,vy,vz          : Eulerian velocity components on regular grid. Must
!                        be of size nd_ set in constructor
!
!-----------------------------------------------------------------
!$  USE threads
    USE mpivars

    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)                      :: this
    INTEGER                                             :: i,j,k,ngp,ngz,nex,nexy,nez
    INTEGER                                             :: nx,nxy,ny,nz
    INTEGER                                             :: jm,km
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: vx,vy,vz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: vxext,vyext,vzext

    ngz  = this%nzghost_
    ngp  = ngz * this%iextperp_
    nexy = (this%nd_(1)+2*ngp) * (this%nd_(2)+2*ngp)
    nex  = this%nd_(1)+2*ngp
    nez  = this%nd_(3)+  ngz
    nx   = this%nd_(1)
    ny   = this%nd_(2)
    nz   = this%nd_(3)
    nxy  = nx*ny

    CALL GPartComm_Copy2Ext(this,vxext,vx)
    CALL GPartComm_Copy2Ext(this,vyext,vy)
    CALL GPartComm_Copy2Ext(this,vzext,vz)
    DO k = 1, ngz  ! bottom extended zones
      km = k-1
!$omp parallel do if(ny.ge.nth) private(jm,i)
      DO j=1,ny
        jm = j-1
        DO i=1,nx
          ! set top bcs:
          vxext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy) = vx(i+(j-1)*nx+(k-1)*nxy)
          vyext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy) = vy(i+(j-1)*nx+(k-1)*nxy)
          vzext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy) = vz(i+(j-1)*nx+(k-1)*nxy)

          ! set bottom bcs:
          vxext(i+ngp+(j+ngp-1)*nex+    (k-1)*nexy) = vx(i+(j-1)*nx+(nz-ngz+k-1)*nxy)
          vyext(i+ngp+(j+ngp-1)*nex+    (k-1)*nexy) = vy(i+(j-1)*nx+(nz-ngz+k-1)*nxy)
          vzext(i+ngp+(j+ngp-1)*nex+    (k-1)*nexy) = vz(i+(j-1)*nx+(nz-ngz+k-1)*nxy)

        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE GPartComm_LocalDataExchMF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_SlabDataExchangeMF(this,vxext,vyext,vzext,vx,vy,vz)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : SlabDataExchangeMF
!  DESCRIPTION: Does bdy exchange of velocity component, vx,vy,vz. Output
!               is to data on extended grids, vxext, vyext, vzexy. 'MF' means
!               that this is the 'multi-field' interface.
!  ARGUMENTS  :
!    this              : 'this' class instance (IN)
!    vxext,vyext,vzext : Eulerian velocity components returned on extended
!                        grid (that used to hold ghost zones)
!    vx,vy,vz          : Eulerian velocity components on regular grid. Must
!                        be of size nd_ set in constructor
!
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)                      :: this
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: vx,vy,vz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: vxext,vyext,vzext

    INTEGER                                             :: itask,j
    IF ( this%intrfc_ .LT. 1 ) THEN
      WRITE(*,*) 'GPartComm_SlabDataExchangeMF: SF interface expected'
      STOP
    ENDIF

    IF ( this%nprocs_ .EQ. 1 ) THEN
      CALL GPartComm_LocalDataExchMF(this,vxext,vyext,vzext,vx,vy,vz)
      RETURN
    ENDIF
    CALL GPartComm_Copy2Ext(this,vxext,vx)
    CALL GPartComm_Copy2Ext(this,vyext,vy)
    CALL GPartComm_Copy2Ext(this,vzext,vz)

    ! Post receives:
    CALL GTStart(this%hcomm_)
    DO j=1,this%nbrcv_  ! from bottom task:
      itask = this%ibrcvp_(j)
      CALL MPI_IRECV(this%rbbuff_(:,j),this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%ibrh_(j),this%ierr_)
    ENDDO
    CALL GTAcc(this%hcomm_)

    CALL GTStart(this%hcomm_)
    DO j=1,this%ntrcv_  ! from top task:
      itask = this%itrcvp_(j)
      CALL MPI_IRECV(this%rtbuff_(:,j),this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%itrh_(j),this%ierr_)
    ENDDO
    CALL GTAcc(this%hcomm_)

    !
    ! send data:
    DO j=1,this%nbsnd_  ! to bottom task:
      itask = this%ibsndp_(j)
      CALL GPartComm_PackMF(this,this%sbbuff_(:,j),vx,vy,vz,j,'b')
      CALL GTStart(this%hcomm_)
      CALL MPI_ISEND(this%sbbuff_(:,j),this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%ibsh_(j),this%ierr_)
      CALL GTAcc(this%hcomm_)
    ENDDO
    DO j=1,this%ntsnd_  ! to top task:
      itask = this%itsndp_(j)
      CALL GPartComm_PackMF(this,this%stbuff_,vx,vy,vz,j,'t')
      CALL GTStart(this%hcomm_)
      CALL MPI_ISEND(this%stbuff_,this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%itsh_(j),this%ierr_)
      CALL GTAcc(this%hcomm_)

    ENDDO

    CALL GTStart(this%hcomm_)
    DO j=1,this%nbsnd_
      CALL MPI_WAIT(this%ibsh_(j),this%istatus_,this%ierr_)
    ENDDO
    DO j=1,this%ntsnd_
      CALL MPI_WAIT(this%itsh_(j),this%istatus_,this%ierr_)
    ENDDO
    DO j=1,this%nbrcv_
      CALL MPI_WAIT(this%ibrh_(j),this%istatus_,this%ierr_)
    ENDDO
    DO j=1,this%ntrcv_
      CALL MPI_WAIT(this%itrh_(j),this%istatus_,this%ierr_)
    ENDDO
    CALL GTAcc(this%hcomm_)

    ! Unpack received data:
    DO j=1,this%nbrcv_
      CALL GPartComm_UnpackMF(this,vxext,vyext,vzext,this%rbbuff_(:,j))
    ENDDO
    DO j=1,this%ntrcv_
      CALL GPartComm_UnpackMF(this,vxext,vyext,vzext,this%rtbuff_(:,j))
    ENDDO

  END SUBROUTINE GPartComm_SlabDataExchangeMF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_PackMF(this,buff,vx,vy,vz,isnd,sdir)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PackMF
!  DESCRIPTION: packs snd buffer with fields; multi-field interface
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    buff    : packed buffer (returned)
!    vx,vy,vz: Eulerian velocity component on regular grid
!              in phys. space (IN)
!    isnd    : which send this is
!    sdir    : 't' for top, 'b' for bottom
!
!-----------------------------------------------------------------
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)          :: this
    INTEGER      ,INTENT   (IN)             :: isnd
    INTEGER                                 :: i,j,k,m,nt,nx,nxy,ny
    INTEGER                                 :: jm,km
    REAL(KIND=GP),INTENT  (OUT),DIMENSION(*):: buff
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*):: vx,vy,vz
    CHARACTER*(*),INTENT   (IN)             :: sdir


    IF ( sdir(1:1).NE.'b' .AND. sdir(1:1).NE.'B' &
    .AND.sdir(1:1).NE.'t' .AND. sdir(1:1).NE.'T' ) THEN
      WRITE(*,*) 'GPartComm_PackMF: Bad direction descriptor'
      STOP
    ENDIF

    nx  = this%nd_(1)
    ny  = this%nd_(2)
    nxy = nx*ny
    IF      ( sdir(1:1) .EQ. 'b' .OR. sdir(1:1) .EQ. 'B' ) THEN
    ! Pack for send to rank at bottom:
    !  ...header
      nt = 1
      buff(1)  = this%ibsndnz_(isnd)       ! no. z-indices included
      DO j = 1, this%ibsndnz_(isnd)
        nt       = nt + 1
        buff(nt) = this%ibsnddst_(isnd,j) ! z-index in extended grid
      ENDDO

    !  ...data
      DO m = 1,this%ibsndnz_(isnd)
        k = this%ibsnd_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            buff(nt) = vx(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

      DO m = 1, this%ibsndnz_(isnd)
        k = this%ibsnd_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            buff(nt) = vy(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

      DO m = 1,this%ibsndnz_(isnd)
        k = this%ibsnd_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            buff(nt) = vz(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

    ELSE !  Pack for send to rank at top:

      ! ...header
      nt = 1
      buff(1)  = this%itsndnz_(isnd)      ! no. z-indices included
      DO j = 1, this%itsndnz_(isnd)
        nt       = nt + 1
        buff(nt) = this%itsnddst_(isnd,j) ! z-index in extended grid
      ENDDO

      ! ...data
      DO m = 1,this%itsndnz_(isnd)
        k = this%itsnd_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            buff(nt) = vx(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

      DO m = 1,this%itsndnz_(isnd)
        k = this%itsnd_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            buff(nt) = vy(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

      DO m = 1,this%itsndnz_(isnd)
        k = this%itsnd_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            buff(nt) = vz(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE GPartComm_PackMF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_UnpackMF(this,vxe,vye,vze,buff)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : UnpackMF
!  DESCRIPTION: Unpacks recv buffer with into extended (single) field
!               Messages received are self-referential, so contain info
!               on where to 'send' recvd data. So, there is no 't' or
!               'b' designation required for unpacking.
!  ARGUMENTS  :
!    this        : 'this' class instance (IN)
!    buff        : packed buffer (input) from which to store into
!                  extended grid quantities.
!    vxe,vye,vze : Eulerian velocity component on extended grid
!                  in phys. space (IN)
!
!-----------------------------------------------------------------
    USE mpivars
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)          :: this
    INTEGER                                 :: i,j,k,m,ngp,ngz,nx,nxy,ny,nz
    INTEGER                                 :: ixy,jm,km,nh
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*):: buff
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*):: vxe,vye,vze

    nx  = this%nd_(1)
    ny  = this%nd_(2)
    nxy = nx*ny
    ngz = this%nzghost_;
    ngp = ngz*this%iextperp_

  ! Unpack from either top or bottom buffer:
    nz = int(buff(1))
    nh = nz + 1 ! no. items in header
    DO m = 1,nz
      k   = int(buff(m+1))
      km  = k-1

      ixy = 1
      DO j = 1, ny
        jm = j-1
        DO i = 1, nx
          vxe(i+ngp+(jm+ngp)*nx+km*nxy) = buff(nh+(m-1)*nxy+ixy)
          ixy = ixy + 1
        ENDDO
      ENDDO

      DO j = 1, ny
        jm = j-1
        DO i = 1, nx
          vye(i+ngp+(jm+ngp)*nx+km*nxy) = buff(nh+(m-1)*nxy+ixy)
          ixy = ixy + 1
        ENDDO
      ENDDO

      DO j = 1, ny
        jm = j-1
        DO i = 1, nx
          vze(i+ngp+(jm+ngp)*nx+km*nxy) = buff(nh+(m-1)*nxy+ixy)
          ixy = ixy + 1
        ENDDO
      ENDDO

    ENDDO

  END SUBROUTINE GPartComm_UnpackMF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_LocalDataExchSF(this,vext,v)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : LocalDataExchSF
!  DESCRIPTION: Does 'bdy exchange' of (single) velocity component, when there's
!               only a single MPI task. This is a single-field interface.
!  ARGUMENTS  :
!    this              : 'this' class instance (IN)
!    vxext,vyext,vzext : Eulerian velocity components returned on extended
!                        grid (that used to hold ghost zones). Only z-conditions
!                        are imposed; lateral periodicity is not handled here.
!    vx,vy,vz          : Eulerian velocity components on regular grid. Must
!                        be of size nd_ set in constructor
!
!-----------------------------------------------------------------
!$  USE threads
    USE mpivars

    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)                      :: this
    INTEGER                                             :: i,j,k,ngp,ngz,nex,nexy,nez
    INTEGER                                             :: nx,nxy,ny,nz
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: v
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: vext

    ngz  = this%nzghost_
    ngp  = ngz * this%iextperp_
    nexy = (this%nd_(1)+2*ngp) * (this%nd_(2)+2*ngp)
    nex  = this%nd_(1)+2*ngp
    nez  = this%nd_(3)+  ngz
    nx   = this%nd_(1)
    ny   = this%nd_(2)
    nz   = this%nd_(3)
    nxy  = nx*ny

    CALL GPartComm_Copy2Ext(this,vext,v)
    DO k = 1, ngz  ! extended zones
!$omp parallel do if(ny.ge.nth) private(i)
      DO j=1,ny
        DO i=1,nx
          ! set top bcs:
          vext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy) = v(i+(j-1)*nx+(k-1)*nxy)

          ! set bottom bcs:
          vext(i+ngp+(j+ngp-1)*nex+    (k-1)*nexy) = v(i+(j-1)*nx+(nz-ngz+k-1)*nxy)

        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE GPartComm_LocalDataExchSF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_SlabDataExchangeSF(this,vext,v)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPartComm_SlabDataExchangeSF
!  DESCRIPTION: Does bdy exchange of field component, v. Output
!               is to data on extended grids, vext. 'SF' means
!               that this is the 'single-field' interface.
!  ARGUMENTS  :
!    this      : 'this' class instance (IN)
!    vext      : Eulerian velocity component returned on extended
!                grid (that used to hold ghost zones in z)
!    v         : Eulerian velocity components on regular grid. Must
!                be of size nd_ set in constructor
!
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)                      :: this
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: v
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: vext

    INTEGER                                             :: itask,j,buffsize

    IF ( this%intrfc_ .GE. 1 ) THEN
      WRITE(*,*) 'GPartComm_SlabDataExchangeSF: MF interface expected'
      STOP
    ENDIF
    IF ( this%nprocs_ .EQ. 1 ) THEN
      CALL GPartComm_LocalDataExchSF(this,vext,v)
      RETURN
    ENDIF

    buffsize = this%nd_(1)*this%nd_(2)*this%nzghost_ + this%nzghost_ + 2

    CALL GPartComm_Copy2Ext(this,vext,v)

    ! post receives:
    CALL GTStart(this%hcomm_)
    DO j=1,this%nbrcv_  ! from bottom task:
      itask = this%ibrcvp_(j)
      CALL MPI_IRECV(this%rbbuff_(:,j),this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%ibrh_(j),this%ierr_)
    ENDDO
    CALL GTAcc(this%hcomm_)


    ! send data:
    DO j=1,this%nbsnd_  ! to bottom task:
      itask = this%ibsndp_(j)
      CALL GPartComm_PackSF(this,this%sbbuff_(:,j),v,j,'b')
      CALL GTStart(this%hcomm_)
      CALL MPI_ISEND(this%sbbuff_(:,j),buffsize,GC_REAL,itask, &
                     1,this%comm_,this%ibsh_(j),this%ierr_)
      CALL GTAcc(this%hcomm_)
    ENDDO
!

    CALL GTStart(this%hcomm_)
    DO j=1,this%ntrcv_  ! from top task:
      itask = this%itrcvp_(j)
      CALL MPI_IRECV(this%rtbuff_(:,j),this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%itrh_(j),this%ierr_)
    ENDDO
    CALL GTAcc(this%hcomm_)


    DO j=1,this%ntsnd_  ! to top task:
      itask = this%itsndp_(j)
      CALL GPartComm_PackSF(this,this%stbuff_(:,j),v,j,'t')
      CALL GTStart(this%hcomm_)
      CALL MPI_ISEND(this%stbuff_(:,j),buffsize,GC_REAL,itask, &
                     1,this%comm_,this%itsh_(j),this%ierr_)
      CALL GTAcc(this%hcomm_)
    ENDDO

    CALL GTStart(this%hcomm_)
    DO j=1,this%nbsnd_
      CALL MPI_WAIT(this%ibsh_(j),this%istatus_,this%ierr_)
    ENDDO
    DO j=1,this%nbrcv_
      CALL MPI_WAIT(this%ibrh_(j),this%istatus_,this%ierr_)
    ENDDO


    DO j=1,this%ntsnd_
      CALL MPI_WAIT(this%itsh_(j),this%istatus_,this%ierr_)
    ENDDO
    DO j=1,this%ntrcv_
      CALL MPI_WAIT(this%itrh_(j),this%istatus_,this%ierr_)
    ENDDO
    CALL GTAcc(this%hcomm_)


    ! Unpack received data:
    DO j=1,this%nbrcv_
      CALL GPartComm_UnpackSF(this,vext,this%rbbuff_(:,j),this%nbuff_,'b',this%ierr_)
    ENDDO
    DO j=1,this%ntrcv_
      CALL GPartComm_UnpackSF(this,vext,this%rtbuff_(:,j),this%nbuff_,'t',this%ierr_)
    ENDDO

  END SUBROUTINE GPartComm_SlabDataExchangeSF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_PackSF(this,buff,v,isnd,sdir)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PPackSF
!  DESCRIPTION: packs snd buffer with (single) field
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    buff    : packed buffer (returned)
!    v       : Eulerian velocity component on regular grid
!              in phys. space (IN)
!    isnd    : which send this is
!    sdir    : 't' for top, 'b' for bottom
!
!-----------------------------------------------------------------
!$  USE threads
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)          :: this
    INTEGER      ,INTENT   (IN)             :: isnd
    INTEGER                                 :: i,j,k,m,nt,nx,ny,nxy
    INTEGER                                 :: jm,km
    REAL(KIND=GP),INTENT  (OUT),DIMENSION(*):: buff
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*):: v
    CHARACTER*(*),INTENT   (IN)             :: sdir


    IF ( sdir(1:1).NE.'b' .AND. sdir(1:1).NE.'B' &
    .AND.sdir(1:1).NE.'t' .AND. sdir(1:1).NE.'T' ) THEN
      WRITE(*,*) 'GPartComm_PackMF: Bad direction descriptor'
      STOP
    ENDIF

    nx  = this%nd_(1)
    ny  = this%nd_(2)
    nxy = nx*ny
    IF      ( sdir(1:1) .EQ. 'b' .OR. sdir(1:1) .EQ. 'B' ) THEN
    ! Pack for send to rank at bottom:
    !  ...header
      nt = 1
      buff(1)  = this%ibsndnz_(isnd)      ! no. z-indices included
      DO j = 1, this%ibsndnz_(isnd)
        nt       = nt + 1
        buff(nt) = this%ibsnddst_(isnd,j) ! z-index in extended grid
      ENDDO

    !  ...data
      DO m = 1,this%ibsndnz_(isnd)
        k = this%ibsnd_(isnd,m)
        km = k-1
!$omp parallel do if(ny.ge.nth) private(jm,i)
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
!            nt = nt + 1
            buff(nt+jm*nx+i) = v(i+jm*nx+km*nxy)
          ENDDO
        ENDDO
        nt = nt + nxy
      ENDDO

    ELSE !  Pack for send to rank at top:

      ! ...header
      nt = 1
      buff(1)  = this%itsndnz_(isnd)      ! no. z-indices included
      DO j = 1, this%itsndnz_(isnd)
        nt       = nt + 1
        buff(nt) = this%itsnddst_(isnd,j) ! z-index in extended grid
      ENDDO

      ! ...data
      DO m = 1,this%itsndnz_(isnd)
        k = this%itsnd_(isnd,m)
        km = k-1
!$omp parallel do if(ny.ge.nth) private(jm,i)
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
!            nt = nt + 1
            buff(nt+jm*nx+i) = v(i+jm*nx+km*nxy)
          ENDDO
        ENDDO
        nt = nt + nxy
      ENDDO

    ENDIF

  END SUBROUTINE GPartComm_PackSF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_UnpackSF(this,vext,buff,nbuff,sb,ierr)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : UnpackSF
!  DESCRIPTION: Unpacks recv buffer with into extended (single) field
!               Messages received are self-referential, so contain info
!               on where to 'send' recvd data. So, there is no 't' or
!               'b' designation required for unpacking.
!  ARGUMENTS  :
!    this        : 'this' class instance (IN)
!    vext        : Eulerian velocity component on extended grid
!                  in phys. space (IN)
!    buff        : packed buffer (input) from which to store into
!                  extended grid quantities.
!    nbuff       : maximum buff size
!    sb          : optional buffer name ,'t' or 'b'.
!    ierr        : err flag: 0 if success; else 1
!
!-----------------------------------------------------------------
!$  USE threads

    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)          :: this
    INTEGER                                 :: i,j,k,m,ngp,nx,nex,nxy,nexy,ny,nz
    INTEGER                                 :: im,ip,ir,ixy,iz,jm,km,nh
    INTEGER      ,INTENT   (IN)             :: nbuff
    INTEGER      ,INTENT(INOUT)             :: ierr ! not used now
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*):: buff
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*):: vext
    CHARACTER(len=1),INTENT(IN),OPTIONAL    :: sb

    ierr = 0
    nx   = this%nd_(1)
    ny   = this%nd_(2)
    nxy  = nx*ny
    ngp  = this%nzghost_*this%iextperp_
    nex  = nx+2*ngp
    nexy = (nx+2*ngp)*(ny+2*ngp)

    ! Unpack from either buffer:
    ! For each task, message is of form:
    !     #z-indices:z-index_0:z-index_1:...:nx*ny_0:nx*ny_1: ...
    nz = int(buff(1))
    nh = nz + 1 ! no. items in header
    DO m = 1, nz
      k   = int(buff(m+1))
      km  = k-1
!      ixy = 1
!$omp parallel do if(ny.ge.nth) private(jm,i,ixy,im,ir)
      DO j = 1, ny
        jm = j-1
        DO i = 1, nx
          ixy = jm*nx + i
          im = i+ngp+(jm+ngp)*nex+km*nexy
          ir = nh+(m-1)*nxy+ixy
          vext(im) = buff(ir)
!          ixy = ixy + 1
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE GPartComm_UnpackSF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_IdentifyTaskV(this,id,pz,nparts,task)
!-----------------------------------------------------------------
!  METHOD     : IdentifyTaskV
!  DESCRIPTION: Identifies particle task owner for scatter. Uses V interface.
!               Note: For this call to work, the particle positions
!               _must_ be periodized on entry.
!
!               This routine is intended to be called  before PartScatterV.
!
!  ARGUMENTS  :
!    this     : 'this' class instance (IN)
!    id       : array of particle ids
!    px,py,px : arrays containing x,y,z positions of particles
!    nparts   : number of particles in pdb
!    newnparts: number of particles in pdb after exchange
!    zmin/max : min/max z-dimensions of current MPI task
!    gext     : (3,2) real array containing global grid extents (start
!               and stop boundaries in each direction).
!-----------------------------------------------------------------
    USE pdbtypes
    USE grid
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)           :: this
    INTEGER      ,INTENT(INOUT)              :: nparts
    INTEGER                                  :: i,j,tsta,tend
    INTEGER      ,INTENT(INOUT),DIMENSION(:) :: id,task
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:) :: pz
    REAL(KIND=GP)                            :: zmaxs(this%nprocs_)
    INTEGER                                  :: nps(this%nprocs_)

    IF ( this%nprocs_ .EQ. 1 ) RETURN ! nothing to do
    
    IF (this%myrank_.EQ.0) THEN
      ! Determine z-boundaries for each task
      DO j = 1,this%nprocs_
        nps(j) = 0
        CALL range(1,nz,this%nprocs_,j-1,tsta,tend)
        zmaxs(j) = REAL(tend,KIND=GP)
      END DO

      ! Determine corresponding task owner for each particle
      DO i = 1,nparts
        j = 1
        DO WHILE (pz(i).GT.zmaxs(j))
          j = j + 1
        END DO
        IF (j.GT.this%nprocs_) THEN
          WRITE(*,*) 'GPartComm_IdentifyTaskV: particle outside of range z =', pz(i)
        END IF
        nps(j) = nps(j) + 1
        task(i) = j - 1
      END DO

      ! Communicate number of particles to send
      ! Build and send particle pack to each task
      DO j = 2,this%nprocs_
        CALL MPI_ISEND(nps(j),1,MPI_INTEGER,j-1,j-1,this%comm_,this%itsh_(1),this%ierr_)
      END DO
    ELSE
      CALL MPI_IRECV(nparts,1,MPI_INTEGER,0,this%myrank_,this%comm_,this%itrh_(1),this%ierr_)
      CALL MPI_WAIT(this%itrh_(1),this%istatus_,this%ierr_)
    END IF

    RETURN

  END SUBROUTINE GPartComm_IdentifyTaskV
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_PartScatterV(this,id,px,py,pz,nparts,task)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PartScatterV
!  DESCRIPTION: Carries out particle scatter from task 0 to the rest.
!               Uses V interface.
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    id      : array of particle ids
!    px,py,px: arrays containing x,y,z positions of particles
!    nparts  : total number of particles in pdb (=maxparts for task 0,
!                                                =0 otherwise)
!    task    : array containing particle owners (tasks)
!-----------------------------------------------------------------
    USE pdbtypes
    USE grid
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)           :: this
    INTEGER      ,INTENT(INOUT)              :: nparts
    INTEGER                                  :: i,j,nsend,nparts0
    INTEGER      ,INTENT(INOUT),DIMENSION(:) :: id
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:) :: px,py,pz
    INTEGER      ,INTENT   (IN),DIMENSION(:) :: task

    IF ( this%nprocs_ .EQ. 1 ) RETURN ! nothing to do

    IF (this%myrank_.EQ.0) THEN
      ! Build and send particle pack to each task
      DO j = 2,this%nprocs_
        nsend = 0
        DO i = 1,this%maxparts_
          IF (task(i).EQ.(j-1)) THEN
            nsend = nsend + 1
            this%itop_(nsend) = i
          END IF
        END DO
        CALL GPartComm_PPackV2(this,this%stbuffp_,id,px,py,pz,this%itop_,nsend)
        CALL MPI_ISEND(this%stbuffp_,nsend,MPI_GPDataPackType,j-1,j-1,this%comm_,this%itsh_(1),this%ierr_)
        CALL MPI_WAIT(this%itsh_(1),this%istatus_,this%ierr_)
      END DO

      ! Concatenate own particles
      nparts0 = 0
      DO i = 1,this%maxparts_
        IF (task(i).EQ.0) THEN
          nparts0 = nparts0 + 1
          id(nparts0) = id(i)
          px(nparts0) = px(i)
          py(nparts0) = py(i)
          pz(nparts0) = pz(i)
        END IF
      END DO
      nparts = nparts0
    ELSE
      CALL MPI_IRECV(this%rtbuffp_,this%maxparts_,MPI_GPDataPackType,0,this%myrank_,this%comm_,this%itrh_(1),this%ierr_)
      CALL MPI_WAIT(this%itrh_(1),this%istatus_,this%ierr_)
      CALL GPartComm_PUnpackV2(this,id,px,py,pz,0,this%rtbuffp_,nparts)
    END IF

    RETURN
 
  END SUBROUTINE GPartComm_PartScatterV
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_PPackV2(this,buff,id,px,py,pz,iind,nind)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PPack
!  DESCRIPTION: Packs send buffer with particles. Uses V interface.
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    buff    : buffer into which to pack particles for sends
!    id      : part. ids
!    px,py,pz: part. locations
!    iind    : pointers into pdb particle arrays for
!              particles to pack
!    nind    : no. particles to pack
!-----------------------------------------------------------------
    USE pdbtypes
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER      ,INTENT   (IN)                :: nind
    INTEGER      ,INTENT   (IN),DIMENSION(:)   :: iind
    INTEGER      ,INTENT   (IN),DIMENSION(:)   :: id
    INTEGER                                    :: j
    TYPE(GPData) ,INTENT(INOUT),DIMENSION(:)   :: buff
    REAL(KIND=GP),INTENT   (IN),DIMENSION(:)   :: px,py,pz

!$omp parallel do if(nind.ge.NMIN_OMP)
    DO j = 1,nind
      buff(j)%idp_   = id(iind(j))
      buff(j)%rp_(1) = px(iind(j))
      buff(j)%rp_(2) = py(iind(j))
      buff(j)%rp_(3) = pz(iind(j))
    ENDDO

  END SUBROUTINE GPartComm_PPackV2
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_PUnpackV2(this,id,px,py,pz,nparts,buff,nbuff)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PUnpackV
!  DESCRIPTION: Unpacks recv buffer with particles. Partlcles
!               will be added directly to the existing particle list.
!               Uses V interface.
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    id      : part. ids
!    px,py,pz: part. locations
!    nparts  : new number of particles in pdb
!              with new particles
!    buff    : buffer from which particle data is read
!    nbuff   : buffer length
!-----------------------------------------------------------------
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER      ,INTENT   (IN)                :: nparts
    INTEGER      ,INTENT   (IN)                :: nbuff
    INTEGER      ,INTENT(INOUT),DIMENSION(:)   :: id
    INTEGER                                    :: j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:)   :: px,py,pz
    TYPE(GPData) ,INTENT   (IN),DIMENSION(:)   :: buff

!$omp parallel do if(nbuff.ge.NMIN_OMP) 
    DO j = 1,nbuff
      id(nparts+j) = buff(j)%idp_
      px(nparts+j) = buff(j)%rp_(1)
      py(nparts+j) = buff(j)%rp_(2)
      pz(nparts+j) = buff(j)%rp_(3)
    ENDDO

  END SUBROUTINE GPartComm_PUnpackV2
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_IdentifyExchV(this,id,pz,nparts,newnparts,zmin,zmax)
!-----------------------------------------------------------------
!  METHOD     : IdentifyExchV
!  DESCRIPTION: Identifies particle Ids for exchange. Uses V interface.
!               Note: For this call to work, the particle positions
!               must _not_ be periodized on entry. In the same way,
!               zmin/zmax must also _not_ be periodized.
!
!               This routine is intended to be called at each stage of
!               an explicit time integration where the particle positions
!               cannot change more than a single zone in x, y, or z
!               in a timestep, before PartExchangeV.
!
!               Note that here, a particle _only on_ zmax is considered
!               to be outside the interval defined by zmin/zmax.
!  ARGUMENTS  :
!    this     : 'this' class instance (IN)
!    id       : array of particle ids
!    px,py,px : arrays containing x,y,z positions of particles
!    nparts   : number of particles in pdb
!    newnparts: number of particles in pdb after exchange
!    zmin/max : min/max z-dimensions of current MPI task
!    gext     : (3,2) real array containing global grid extents (start
!               and stop boundaries in each direction).
!-----------------------------------------------------------------
    USE pdbtypes
    USE grid
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)           :: this
    INTEGER      ,INTENT(INOUT)              :: nparts,newnparts
    INTEGER                                  :: j,ibrank,itrank
    INTEGER                                  :: nrtop,nrbot
    INTEGER      ,INTENT(INOUT),DIMENSION(:) :: id
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:) :: pz
    REAL(KIND=GP),INTENT   (IN)              :: zmin,zmax

    IF ( this%nprocs_ .EQ. 1 ) RETURN ! nothing to do

    itrank = modulo(this%myrank_+1,this%nprocs_)
    ibrank = this%myrank_-1
    IF ( ibrank.LT.0 ) ibrank = this%nprocs_-1

    ! Find pointers into particle lists for parts that
    ! must be sent to the top and bottom tasks:
    this%nbot_ = 0
    this%ntop_ = 0
    DO j = 1, nparts
      IF ( pz(j).LT.zmin ) THEN ! bottom
        this%nbot_ = this%nbot_ + 1
        this%ibot_(this%nbot_) = j
      ELSE IF ( pz(j).GE.zmax ) THEN ! top
        this%ntop_ = this%ntop_ + 1
        this%itop_(this%ntop_) = j
      ENDIF
    ENDDO

    ! Post receives:
    CALL GTStart(this%hcomm_)
    CALL MPI_IRECV(nrbot,1,MPI_INTEGER,ibrank,1,this%comm_,this%ibrh_(1),this%ierr_)
    CALL MPI_IRECV(nrtop,1,MPI_INTEGER,itrank,1,this%comm_,this%itrh_(1),this%ierr_)
    CALL GTAcc(this%hcomm_)

    ! send data:
    CALL MPI_ISEND(this%nbot_,1,MPI_INTEGER,ibrank,1,this%comm_,this%ibsh_(1),this%ierr_)
    CALL MPI_ISEND(this%ntop_,1,MPI_INTEGER,itrank,1,this%comm_,this%itsh_(1),this%ierr_)
 
    CALL MPI_WAIT(this%ibrh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%itrh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%ibsh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%itsh_(1),this%istatus_,this%ierr_)
    CALL GTAcc(this%hcomm_)

    newnparts = nparts - (this%ntop_+this%nbot_) + (nrbot+nrtop)
 
    RETURN

  END SUBROUTINE GPartComm_IdentifyExchV
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_PartExchangeV(this,id,px,py,pz,nparts,zmin,zmax,stg)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PartExchangeV
!  DESCRIPTION: Carries out particle exchange. Particles will
!               be re-ordered after this call. Uses V interface.
!               Must be called _after_ IdentifyExchV.
!               Note: For this call to work, the particle positions
!               must _not_ be periodized on entry. In the same way,
!               zmin/zmax must also _not_ be periodized.
!
!               This routine is intended to be called at each stage of
!               an explicit time integration where the particle positions
!               cannot change more than a single zone in x, y, or z
!               in a timestep.
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    id      : array of particle ids
!    px,py,px: arrays containing x,y,z positions of particles
!    nparts  : number of particles in pdb
!    zmin/max: min/max z-dimensions of current MPI task
!    gext    : (3,2) real array containing global grid extents (start and
!              stop boundaries in each direction).
!    stg     : communicator stage (INIT, UPDT, END, UNIQ)
!-----------------------------------------------------------------
    USE pdbtypes
    USE grid
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)           :: this
    INTEGER      ,INTENT(INOUT)              :: nparts
    INTEGER                                  :: j,ibrank,itrank,ng,nrt,nrb
    INTEGER      ,INTENT(INOUT),DIMENSION(:) :: id
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:) :: px,py,pz
    REAL(KIND=GP),INTENT   (IN)              :: zmin,zmax
    INTEGER      ,INTENT   (IN), OPTIONAL    :: stg
    INTEGER                                  :: stage

    IF ( this%nprocs_ .EQ. 1 ) RETURN ! nothing to do

    itrank = modulo(this%myrank_+1,this%nprocs_)
    ibrank = this%myrank_-1
    IF ( ibrank.LT.0 ) ibrank = this%nprocs_-1

    IF (PRESENT(stg)) THEN
       stage = stg
    ELSE
       stage = GPEXCH_UNIQ
    ENDIF

    IF ((stage.NE.GPEXCH_END).AND.(stage.NE.GPEXCH_UNIQ)) THEN
       ng = nparts
!$omp parallel do if(ng.ge.NMIN_OMP)
       DO j = 1,ng
         this%oldid_(j) = id(j)
       END DO
    ENDIF

    ! Post receives:
    CALL GTStart(this%hcomm_)
!    CALL MPI_IRECV(this%rbbuff_(:,1),2*this%nbuff_,GC_REAL,ibrank, &
!                   1,this%comm_,this%ibrh_(1),this%ierr_)
!    CALL MPI_IRECV(this%rtbuff_(:,1),2*this%nbuff_,GC_REAL,itrank, &
!                   1,this%comm_,this%itrh_(1),this%ierr_)
    CALL GTAcc(this%hcomm_)

    !
    ! send data:
!    CALL GPartComm_PPackV(this,this%sbbuff_(:,1),this%nbuff_,id,px,py,pz,nparts,this%ibot_,this%nbot_)
    CALL GPartComm_PPackV2(this,this%sbbuffp_,id,px,py,pz,this%ibot_,this%nbot_)
    CALL GTStart(this%hcomm_)
!    CALL MPI_ISEND(this%sbbuff_,4*this%nbot_+1,GC_REAL,ibrank, &
!                   1,this%comm_,this%ibsh_(1),this%ierr_)
    CALL MPI_ISEND(this%sbbuffp_,this%nbot_,MPI_GPDataPackType,ibrank, &
                   1,this%comm_,this%ibsh_(1),this%ierr_)
    CALL GTAcc(this%hcomm_)

!    CALL GPartComm_PPackV(this,this%stbuff_(:,1),this%nbuff_,id,px,py,pz,nparts,this%itop_,this%ntop_)
    CALL GPartComm_PPackV2(this,this%stbuffp_,id,px,py,pz,this%itop_,this%ntop_)
    CALL GTStart(this%hcomm_)
!    CALL MPI_ISEND(this%stbuff_,4*this%ntop_+1,GC_REAL,itrank, &
!                   1,this%comm_,this%itsh_(1),this%ierr_)
    CALL MPI_ISEND(this%stbuffp_,this%ntop_,MPI_GPDataPackType,itrank, &
                   1,this%comm_,this%itsh_(1),this%ierr_)
    CALL GTAcc(this%hcomm_)

    ! Concatenate partcle list to remove particles sent away:
    CALL GPartComm_ConcatV(this,id,px,py,pz,nparts,this%ibot_,this%nbot_,this%itop_,this%ntop_)

    CALL GTStart(this%hcomm_)
!    CALL MPI_WAIT(this%ibrh_(1),this%istatus_,this%ierr_)
!    CALL MPI_WAIT(this%itrh_(1),this%istatus_,this%ierr_)
    CALL MPI_RECV(this%rbbuffp_,this%maxparts_,MPI_GPDataPackType,ibrank, &
                  1,this%comm_,this%ibrh_(1),this%ierr_)
    CALL MPI_GET_COUNT(this%ibrh_,MPI_GPDataPackType,nrb)
    CALL MPI_RECV(this%rtbuffp_,this%maxparts_,MPI_GPDataPackType,itrank, &
                  1,this%comm_,this%itrh_(1),this%ierr_)
    CALL MPI_GET_COUNT(this%itrh_,MPI_GPDataPackType,nrt)
    CALL MPI_WAIT(this%ibsh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%itsh_(1),this%istatus_,this%ierr_)
    CALL GTAcc(this%hcomm_)

    ! Update particle list:
!    CALL GPartComm_PUnpackV(this,id,px,py,pz,nparts,this%rbbuff_(:,1),this%nbuff_)
!    CALL GPartComm_PUnpackV(this,id,px,py,pz,nparts,this%rtbuff_(:,1),this%nbuff_)
    CALL GPartComm_PUnpackV2(this,id,px,py,pz,nparts,this%rbbuffp_,nrb)
    nparts = nparts + nrb
    CALL GPartComm_PUnpackV2(this,id,px,py,pz,nparts,this%rtbuffp_,nrt)
    nparts = nparts + nrt
    
    IF ((stage.NE.GPEXCH_END).AND.(stage.NE.GPEXCH_UNIQ)) THEN
!$omp parallel do if(ng.ge.NMIN_OMP)
       DO j = 1,ng
         id(j) = this%oldid_(j)
       END DO
       nparts = ng
    ENDIF

  END SUBROUTINE GPartComm_PartExchangeV
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_PPackV(this,buff,nbuff,id,px,py,pz,nparts,iind,nind)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PPack
!  DESCRIPTION: Packs send buffer with particles. Uses V interface.
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    buff    : buffer into which to pack particles for sends
!    nbuff   : max buffer length
!    id      : part. ids
!    px,py,pz: part. locations
!    nparts  : number of particles in pdb
!    iind    : pointers into pdb particle arrays for
!              particles to pack
!    nind    : no. particles to pack
!-----------------------------------------------------------------
    USE pdbtypes
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER      ,INTENT   (IN)                :: nbuff,nparts,nind
    INTEGER      ,INTENT   (IN),DIMENSION(:)   :: iind
    INTEGER      ,INTENT   (IN),DIMENSION(:)   :: id
    INTEGER                                    :: j,nb
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:)   :: buff
    REAL(KIND=GP),INTENT   (IN),DIMENSION(:)   :: px,py,pz

    buff(1) = nind
!    nb = 1
!$omp parallel do if(nind.ge.NMIN_OMP) private(nb)
    DO j = 1, nind
      nb = 4*j - 3
      buff(nb+1) = id(iind(j))
      buff(nb+2) = px(iind(j))
      buff(nb+3) = py(iind(j))
      buff(nb+4) = pz(iind(j))
!      nb = nb + 4
    ENDDO

  END SUBROUTINE GPartComm_PPackV
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_ConcatV(this,id,px,py,pz,nparts,ibind,nbind,itind,ntind)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : ConcatV
!  DESCRIPTION: Removes particles at indices itind,ibind,and
!               concatenates the particles list, using V interface
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    buff    : buffer into which to pack particles for sends
!    id      : part. ids
!    px,py,pz: part. locations
!    nparts  : number of particles into pdb
!              updated
!    ibind   : list of indices of parts sent to bottom task
!    nbind   : no. indices in ibind
!    itind   : list of indices of parts sent to top
!    ntind   : no. indices in itind
!-----------------------------------------------------------------
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)                :: this
    INTEGER      ,INTENT(INOUT)                   :: nparts
    INTEGER      ,INTENT   (IN)                   :: nbind,ntind
    INTEGER      ,INTENT   (IN),DIMENSION(:)      :: ibind,itind
    INTEGER      ,INTENT(INOUT),DIMENSION(:)      :: id
    INTEGER                                       :: i,j,ngood
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:)      :: px,py,pz

    IF ((nbind+ntind).EQ.0) RETURN ! Nothing to do if no particle left

    DO j = 1, nbind
      id(ibind(j)) = GPNULL
    ENDDO
    DO j = 1, ntind
      id(itind(j)) = GPNULL
    ENDDO

    ngood = nparts - (nbind+ntind)
    j     = 1
    DO i = 1, ngood
      DO WHILE (id(j).EQ.GPNULL)
         j = j + 1
      ENDDO
      id(i) = id(j)
      px(i) = px(j)
      py(i) = py(j)
      pz(i) = pz(j)
      j = j + 1
    ENDDO
    nparts = ngood

  END SUBROUTINE GPartComm_ConcatV
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_PUnpackV(this,id,px,py,pz,nparts,buff,nbuff)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PUnpackV
!  DESCRIPTION: Unpacks recv buffer with particles. Partlcles
!               will be added directly to the existing particle list.
!               Uses V interface.
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    id      : part. ids
!    px,py,pz: part. locations
!    nparts  : new number of particles in pdb
!              with new particles
!    buff    : buffer from which particle data is read
!    nbuff   : buffer length
!-----------------------------------------------------------------
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)                :: this
    INTEGER      ,INTENT(INOUT)                   :: nparts
    INTEGER      ,INTENT   (IN)                   :: nbuff
    INTEGER      ,INTENT(INOUT),DIMENSION(:)      :: id
    INTEGER                                       :: j,nb,nind
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:)      :: px,py,pz
    REAL(KIND=GP),INTENT   (IN),DIMENSION(:)      :: buff

!    nb = 1
     nind = int(buff(1))
!$omp parallel do if(nind.ge.NMIN_OMP) private(nb)
    DO j = 1,nind
!      nparts = nparts + 1
      nb = 4*j - 3
      id(nparts+j) = int(buff(nb+1))
      px(nparts+j) =     buff(nb+2)
      py(nparts+j) =     buff(nb+3)
      pz(nparts+j) =     buff(nb+4)
!      nb = nb+4
    ENDDO
    nparts = nparts + nind

  END SUBROUTINE GPartComm_PUnpackV
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_PartExchangePDB(this,pdb,nparts,zmin,zmax)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PartExchangePDB
!  DESCRIPTION: Carries out particle exchange. Particles will
!               be re-ordered after this call. Uses PDB interface.
!               Note: For this call to work, the particle positions
!               must be periodized on entry. In the same way,
!               zmin/zmax must also be periodized.
!
!               This routine is intended to be called at each stage of
!               an explicit time integration where the particle positions
!               cannot change more than a single zone in x, y, or z
!               in a timestep.
!
!               Note that here, a particle _on_ either zmin or zmax
!               is considered to be outside the interval defined
!               by zmin/zmax.
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    pdb     : part. d.b.
!    nparts  : number of particles in pdb
!    zmin/max: min/max z-dimensions of current MPI task
!-----------------------------------------------------------------
    USE pdbtypes
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER      ,INTENT(INOUT)                :: nparts
    INTEGER                                    :: j,ibrank,itrank
    TYPE(GPDBrec),INTENT(INOUT),DIMENSION(*)   :: pdb
    REAL(KIND=GP),INTENT   (IN)                :: zmin,zmax

    IF ( this%nprocs_ .EQ. 1 ) RETURN ! nothing to do

    itrank = modulo(this%myrank_,this%nprocs_)
    ibrank = this%myrank_-1
    IF ( ibrank.LT.0 ) ibrank = this%nprocs_-1

    ! Find pointers into particle lists for parts that must
    ! be sent to the top and bottom tasks:
    this%nbot_ = 0
    this%ntop_ = 0
    IF ( this%myrank_ .EQ. 0 ) THEN

      DO j = 0, nparts
        IF ( pdb(j)%z.GE.zmax ) THEN ! bottom
          this%nbot_ = this%nbot_ + 1
          this%ibot_(this%nbot_) = j
        ELSE
          this%ntop_ = this%ntop_ + 1
          this%itop_(this%ntop_) = j
        ENDIF
      ENDDO

    ELSE IF ( this%myrank_ .EQ. this%nprocs_-1) THEN

      DO j = 0, nparts
        IF ( pdb(j)%z.LE.zmin ) THEN ! top
          this%ntop_ = this%ntop_ + 1
          this%itop_(this%ntop_) = j
        ELSE
          this%nbot_ = this%nbot_ + 1
          this%ibot_(this%nbot_) = j
        ENDIF
      ENDDO

    ELSE

    DO j = 0, nparts
        IF ( pdb(j)%z.LE.zmin ) THEN ! bottom
        this%nbot_ = this%nbot_ + 1
          this%ibot_(this%nbot_) = j
        ENDIF
        IF ( pdb(j)%z.GE.zmax ) THEN ! top
          this%ntop_ = this%ntop_ + 1
          this%itop_(this%ntop_) = j
        ENDIF
      ENDDO

    ENDIF

    ! Post receives:
    CALL GTStart(this%hcomm_)
    CALL MPI_IRECV(this%rbbuff_,this%nbuff_,GC_REAL,ibrank, &
                   1,this%comm_,this%ibrh_(1),this%ierr_)
    CALL MPI_IRECV(this%rtbuff_,this%nbuff_,GC_REAL,itrank, &
                   1,this%comm_,this%itrh_(1),this%ierr_)
    CALL GTAcc(this%hcomm_)


    !
    ! send data:
    CALL GPartComm_PPackPDB(this,this%sbbuff_,this%nbuff_,pdb,nparts,this%ibot_,this%nbot_)
    CALL GTStart(this%hcomm_)
    CALL MPI_ISEND(this%sbbuff_,this%nbuff_,GC_REAL,ibrank, &
                   1,this%comm_,this%itsh_(1),this%ierr_)
    CALL GTAcc(this%hcomm_)

    CALL GPartComm_PPackPDB(this,this%sbbuff_,this%nbuff_,pdb,nparts,this%itop_,this%ntop_)
    CALL GTStart(this%hcomm_)
    CALL MPI_ISEND(this%stbuff_,this%nbuff_,GC_REAL,itrank, &
                   1,this%comm_,this%itsh_(1),this%ierr_)
    CALL GTAcc(this%hcomm_)


    ! Concatenate partcle list to remove particles sent away:
    CALL GPartComm_ConcatPDB(this,pdb,nparts,this%ibot_,&
                          this%nbot_,this%itop_,this%ntop_)

    CALL GTStart(this%hcomm_)
    CALL MPI_WAIT(this%ibrh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%ibrh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%ibsh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%itsh_(1),this%istatus_,this%ierr_)
    CALL GTAcc(this%hcomm_)


    ! Update particle list:
    CALL GPartComm_PUnpackPDB(this,pdb,nparts,this%rbbuff_,this%nbuff_)
    CALL GPartComm_PUnpackPDB(this,pdb,nparts,this%rtbuff_,this%nbuff_)

  END SUBROUTINE GPartComm_PartExchangePDB
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_PPackPDB(this,buff,nbuff,pdb,nparts,iind,nind)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PPackPDB
!  DESCRIPTION: Packs send buffer with particles. Uses PDB interface.
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    buff    : buffer into which to pack particles for sends
!    nbuff   : max buffer length
!    pdb     : part. d.b.
!    nparts  : number of particles in pdb
!    iind    : pointers into pdb particle arrays for
!              particles to pack
!    nind    : no. particles to pack
!-----------------------------------------------------------------
    USE pdbtypes
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER      ,INTENT(INOUT)                :: nbuff,nparts,nind
    INTEGER      ,INTENT(INOUT),DIMENSION(*)   :: iind
    INTEGER                                    :: j,nb
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)   :: buff
    TYPE(GPDBrec),INTENT(INOUT),DIMENSION(*)   :: pdb

    buff(1) = nind
    nb = 1
    DO j = 1, nind
      buff(nb+1) = pdb(iind(j))%id
      buff(nb+2) = pdb(iind(j))%x
      buff(nb+3) = pdb(iind(j))%y
      buff(nb+4) = pdb(iind(j))%z
      nb = nb + 4
    ENDDO

  END SUBROUTINE GPartComm_PPackPDB
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_ConcatPDB(this,pdb,nparts,ibind,nbind,itind,ntind)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : ConcatPDB
!  DESCRIPTION: Removes particles at indices itind,ibind,and
!               concatenates the particles list, using PDB interface
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    buff    : buffer into which to pack particles for sends
!    pdb     : part. d.b.
!    nparts  : number of particles into pdb
!              updated
!    ibind   : list of indices of parts sent to bottom task
!    nbind   : no. indices in ibind
!    itind   : list of indices of parts sent to top
!    ntind   : no. indices in itind
!-----------------------------------------------------------------
    USE pdbtypes
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER      ,INTENT(INOUT)                :: nparts
    INTEGER      ,INTENT   (IN)                :: nbind,ntind
    INTEGER      ,INTENT   (IN),DIMENSION(*)   :: ibind,itind
    INTEGER                                    :: i,j,ngood
    TYPE(GPDBrec),INTENT(INOUT),DIMENSION(*)   :: pdb

    DO j = 1, nbind
      pdb(ibind(j))%id = GPNULL
    ENDDO
    DO j = 1, nbind
      pdb(itind(j))%id = GPNULL
    ENDDO

    ngood = nparts - (nbind+ntind)
    j     = 1
    DO i = 1, ngood
      DO WHILE ( j.LE.nparts .AND. pdb(j)%id.EQ.GPNULL )
        j = j + 1
      ENDDO
      IF ( j.LE.nparts .AND. j.NE.i ) THEN
        pdb(i)%id = pdb(j)%id; pdb(j)%id = GPNULL
        pdb(i)%x = pdb(j)%x
        pdb(i)%y = pdb(j)%y
        pdb(i)%z = pdb(j)%z
      ENDIF

    ENDDO
    nparts = ngood

  END SUBROUTINE GPartComm_ConcatPDB
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_PUnpackPDB(this,pdb,nparts,buff,nbuff)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PUnpackPDB
!  DESCRIPTION: Unpacks recv buffer with particles. Partlcles
!               will be added directly to the existing particle list.
!               Uses PDB interface.
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    pdb     : part. d.b.
!    nparts  : new number of particles in pdb
!              with new particles
!    buff    : buffer from which particle data is read
!    nbuff   : buffer length
!-----------------------------------------------------------------
    USE pdbtypes
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER      ,INTENT(INOUT)                :: nparts
    INTEGER      ,INTENT   (IN)                :: nbuff
    INTEGER                                    :: j,nb
    TYPE(GPDBrec),INTENT(INOUT),DIMENSION(*)   :: pdb
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)   :: buff

    nb = 1
    DO j = 1, int(buff(1))
      nparts = nparts + 1
      pdb(nparts)%id = int(buff(nb+1))
      pdb(nparts)%x =      buff(nb+2)
      pdb(nparts)%y =      buff(nb+3)
      pdb(nparts)%z =      buff(nb+4)
      nb = nb+4
    ENDDO

  END SUBROUTINE GPartComm_PUnpackPDB
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_Transpose(this,ofield,od,ifield,id,rank,tmp)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Transpose
!  DESCRIPTION: Does global transpose to take a x-y complete field,
!               infield, to a yz-complete field, outfield. Handles
!               2D and 3D fields.
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    ofield  : output field, yz-complete
!    od      : local dimensions of ofield.
!    ifield  : input field that is xy complete
!    id      : local dimensions of ifield.
!    rank    : rank of field (how many 'od, id' array elements)
!    tmp     : real field of size required to hold field
!              transpose locally (i.e., of size ofield)
!-----------------------------------------------------------------
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER      ,INTENT(INOUT),DIMENSION(3,2) :: od,id
    INTEGER      ,INTENT   (IN)                :: rank
    INTEGER                                    :: i,ii,j,jj,k,kk
    INTEGER                                    :: igetfrom,iproc,irank,isendto,istrip
    REAL(KIND=GP),INTENT(INOUT)                :: &
      ofield(od(1,1):od(1,2),od(2,1):od(2,2),od(3,1):od(3,2))
    REAL(KIND=GP),INTENT(INOUT)                :: &
      ifield(id(1,1):id(1,2),id(2,1):id(2,2),id(3,1):id(3,2))
    REAL(KIND=GP),INTENT(INOUT)                :: &
      tmp   (od(3,1):od(3,2),od(2,1):od(2,2),od(1,1):od(1,2))

    IF ( .NOT.this%btransinit_ ) THEN
      IF ( rank.EQ.2 ) THEN
        CALL GPartComm_InitTrans2D(this)
      ENDIF
      IF ( rank.EQ.3 ) THEN
        CALL GPartComm_InitTrans3D(this)
      ENDIF
    ENDIF

    ! NOTE: rank is transpose problem rank; irank is MPI rank...

    CALL GTStart(this%hcomm_)
    DO iproc = 0, this%nprocs_-1, this%nstrip_
       DO istrip=0, this%nstrip_-1
          irank = iproc + istrip

          isendto = this%myrank_ + irank
          IF ( isendto .ge. this%nprocs_ ) isendto = isendto - this%nprocs_

          igetfrom = this%myrank_- irank
          IF ( igetfrom .lt. 0 ) igetfrom = igetfrom + this%nprocs_
          CALL MPI_IRECV(tmp,1,this%itypeip_(igetfrom),igetfrom,      &
                        1,this%comm_,this%igrh_(irank),this%ierr_)

          IF ( this%ierr_ .ne. mpi_success ) THEN
            WRITE(*,*)'Transpose: irecv ierr=',this%ierr_
            STOP
          endif
          CALL MPI_ISEND(ifield,1,this%itypekp_(isendto),isendto, &
                        1,this%comm_,this%igsh_(irank),this%ierr_)
          IF ( this%ierr_ .ne. mpi_success ) THEN
            WRITE(*,*)'Transpose: isnd ierr=',this%ierr_
            STOP
          ENDIF
       ENDDO

       DO istrip=0, this%nstrip_-1
          irank = iproc + istrip
          CALL MPI_WAIT(this%igsh_(irank),this%istatus_,this%ierr_)
          IF ( this%ierr_ .ne. mpi_success ) THEN
            WRITE(*,*)'Transpose: Send Wait: ierr=',this%ierr_
            STOP
          ENDIF
          CALL MPI_WAIT(this%igrh_(irank),this%istatus_,this%ierr_)
          IF ( this%ierr_ .ne. mpi_success ) THEN
            WRITE(*,*)'Transpose: Rcv Wait: ierr=',this%ierr_
            STOP
          ENDIF
       ENDDO
    ENDDO
    CALL GTAcc(this%hcomm_)

    IF ( rank .EQ. 3 ) THEN

!!!$omp parallel do if ((idims(3)-1)/this%csize_.ge.this%nth_) private (jj,kk,i,j,k)
     DO ii = od(3,1),od(3,2),this%csize_
!!!$omp parallel do if ((idims(3)-1)/this%csize_.lt.this%nth_) private (kk,i,j,k)
        DO jj = od(2,1),od(2,2),this%csize_
           DO kk = od(1,1),od(1,2),this%csize_

              DO i = ii,min(od(3,2)-od(3,1)+1,ii+this%csize_-1)
                DO j = jj,min(od(2,2)-od(2,1)+1,jj+this%csize_-1)
                  DO k = kk,min(od(1,2)-od(1,1)+1,kk+this%csize_-1)
                     ofield(k,j,i) = tmp(i,j,k)
                  END DO
                END DO
              END DO

           END DO
        END DO
     END DO

    ELSE

      write(*,*) 'GPartComm_Transpose: rank two not implemented'
      stop
    ENDIF

  END SUBROUTINE GPartComm_Transpose
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_ITranspose(this,ofield,od,ifield,id,rank,tmp)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : ITranspose
!  DESCRIPTION: Does global 'inverse'transpose to take a x-y complete field,
!               infield, to a yz-complete field, outfield. Handles
!               2D and 3D fields.
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    ofield  : output field, yz-complete 
!    od      : local dimensions of ofield. 
!    ifield  : input field that is xy complete
!    id      : local dimensions of ifield. 
!    rank    : rank of field (how many 'od, id' array elements)
!    tmp     : real field of size required to hold field locally
!              (i.e., of size ifield)
!-----------------------------------------------------------------
    USE gtimer
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER      ,INTENT(INOUT),DIMENSION(3,2) :: od,id
    INTEGER      ,INTENT   (IN)                :: rank
    INTEGER                                    :: i,ii,j,jj,k,kk
    INTEGER                                    :: igetfrom,iproc,irank,isendto,istrip
    INTEGER                                    :: nx,ny,nz,nxy,nzy
    REAL(KIND=GP),INTENT(INOUT)                :: &
      ofield(od(1,1):od(1,2),od(2,1):od(2,2),od(3,1):od(3,2))
    REAL(KIND=GP),INTENT(INOUT)                :: &
      ifield(id(1,1):id(1,2),id(2,1):id(2,2),id(3,1):id(3,2))
    REAL(KIND=GP),INTENT(INOUT)                :: &
      tmp   (id(3,1):id(3,2),id(2,1):id(2,2),id(1,1):id(1,2))

    IF ( .NOT.this%btransinit_ ) THEN
      IF ( rank.EQ.2 ) THEN
        CALL GPartComm_InitTrans2D(this)
      ENDIF
      IF ( rank.EQ.3 ) THEN
        CALL GPartComm_InitTrans3D(this)
      ENDIF
    ENDIF

    ! NOTE: rank is transpose problem rank; irank is MPI rank...

    IF ( rank .EQ. 3 ) THEN

!!!$omp parallel do if ((idims(3)-1)/this%csize_.ge.this%nth_) private (jj,kk,i,j,k)
     DO ii = id(3,1),id(3,2),this%csize_
!!!$omp parallel do if ((idims(3)-1)/this%csize_.lt.this%nth_) private (kk,i,j,k)
        DO jj = id(2,1),id(2,2),this%csize_
           DO kk = id(1,1),id(1,2),this%csize_

              DO i = ii,min(id(3,2)-id(3,1)+1,ii+this%csize_-1)
                DO j = jj,min(id(2,2)-id(2,1)+1,jj+this%csize_-1)
                  DO k = kk,min(id(1,2)-id(1,1)+1,kk+this%csize_-1)
                     tmp(i,j,k) = ifield(k,j,i)
                  END DO
                END DO
              END DO

           END DO
        END DO
     END DO

    ELSE

    ENDIF

    CALL GTStart(this%hcomm_)
    DO iproc = 0, this%nprocs_-1, this%nstrip_
       DO istrip=0, this%nstrip_-1
          irank = iproc + istrip

          isendto = this%myrank_ + irank
          IF ( isendto .ge. this%nprocs_ ) isendto = isendto - this%nprocs_

          igetfrom = this%myrank_- irank
          IF ( igetfrom .lt. 0 ) igetfrom = igetfrom + this%nprocs_
          CALL MPI_IRECV(ofield,1,this%itypekp_(igetfrom),igetfrom,      &
                        1,this%comm_,this%igrh_(irank),this%ierr_)

          IF  ( this%ierr_ .ne. mpi_success ) THEN
            WRITE (*,*)'ITranspose: irecv ierr=',this%ierr_
            STOP
         ENDIF
          CALL MPI_ISEND(tmp,1,this%itypeip_(isendto),isendto, &
                        1,this%comm_,this%igsh_(irank),this%ierr_)
          IF ( this%ierr_ .ne. mpi_success ) THEN
            WRITE (*,*)'ITranspose: isnd ierr=',this%ierr_
            STOP
          ENDIF
       ENDDO

       DO istrip=0, this%nstrip_-1
          irank = iproc + istrip
          CALL MPI_WAIT(this%igsh_(irank),this%istatus_,this%ierr_)
          IF ( this%ierr_ .ne. mpi_success ) THEN
            WRITE(*,*)'ITranspose: Send Wait: ierr=',this%ierr_
            STOP
          endif
          CALL MPI_WAIT(this%igrh_(irank),this%istatus_,this%ierr_)
          IF ( this%ierr_ .ne. mpi_success ) THEN
            WRITE (*,*)'ITranspose: Rcv Wait: ierr=',this%ierr_
            STOP
          ENDIF
       ENDDO
    ENDDO
    CALL GTAcc(this%hcomm_)

  END SUBROUTINE GPartComm_ITranspose
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_InitTrans2D(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : InitTranspose2D
!  DESCRIPTION: Initializes communcation quantities for 2D transpose.
!               Derived from 2D/src/fftp-3/fftp2d.fpp:fftp2d_create_block
!               and calls function from that module.
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!-----------------------------------------------------------------
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER                                    :: ista,iend
    INTEGER                                    :: jsta,jend
    INTEGER                                    :: irank,jrank
    INTEGER                                    :: itemp1,itemp2

    write(*,*)'GPartComm_InitTrans2D: block2d not resolved'
    stop

    CALL range(1,this%nd_(2),this%nprocs_,this%myrank_,jsta,jend)
    DO irank = 0,this%nprocs_-1
       CALL range(1,this%nd_(1),this%nprocs_,irank,ista,iend)
!      CALL block2d(1,this%nd_(1),jsta,ista,iend,jsta,jend, &
!                   GC_REAL,itemp1)
       this%itypekp_(irank) = itemp1
    END DO
    CALL range(1,this%nd_(1),this%nprocs_,this%myrank_,ista,iend)
    DO jrank = 0,this%nprocs_-1
       CALL range(1,this%nd_(2),this%nprocs_,jrank,jsta,jend)
!      CALL block2d(ista,iend,1,ista,iend,jsta,jend,  &
!                  GC_REAL,itemp2)
       this%itypeip_(jrank) = itemp2
    END DO
    this%btransinit_ = .TRUE.

    RETURN

  END SUBROUTINE GPartComm_InitTrans2D
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_InitTrans3D(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : InitTranspose3D
!  DESCRIPTION: Initializes communcation quantities for 3D transpose.
!               Derived from 3D/src/fftp-3/fftp3d.fpp:fftp3d_create_block
!               and calls function from that module.
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!-----------------------------------------------------------------
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER                                    :: ista,iend
    INTEGER                                    :: ksta,kend
    INTEGER                                    :: irank,krank
    INTEGER                                    :: itemp1,itemp2


    CALL range(1,this%nd_(3),this%nprocs_,this%myrank_,ksta,kend)
    DO irank = 0,this%nprocs_-1
       CALL range(1,this%nd_(1),this%nprocs_,irank,ista,iend)
       CALL block3d(1,this%nd_(1),1,this%nd_(2),ksta,ista,iend, &
                    1,this%nd_(2),ksta,kend,GC_REAL,itemp1)
       this%itypekp_(irank) = itemp1
    END DO
    CALL range(1,this%nd_(1),this%nprocs_,this%myrank_,ista,iend)
    DO krank = 0,this%nprocs_-1
       CALL range(1,this%nd_(3),this%nprocs_,krank,ksta,kend)
       CALL block3d(ista,iend,1,this%nd_(2),1,ista,iend, &
                   1,this%nd_(2),ksta,kend,GC_REAL,itemp2)
       this%itypeip_(krank) = itemp2
    END DO
    this%btransinit_ = .TRUE.

    RETURN

  END SUBROUTINE GPartComm_InitTrans3D
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  FUNCTION GPartComm_GetNumGhost(this) result(nzghost_result)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetNumGhost
!  DESCRIPTION: Get no. ghost zones expected to be transferred.
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!-----------------------------------------------------------------
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER                                    :: nzghost_result

    nzghost_result = this%nzghost_

  END FUNCTION GPartComm_GetNumGhost
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_VDBSynch_t0(this,gvdb,ngvdb,id,lx,ly,lz,nl)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : VDBSynch_t0
!  DESCRIPTION: Synch up only task 0 VDB from local vector data
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    gvdb    : VDB containing part. position records, returned (task 0 only)
!    ngvdb   : no. records in global VDB. Fixed on entry
!    id      : local part. ids
!    lx,ly,lz: local part. d.b. vectors
!    nl      : no. parts. in local pdb
!-----------------------------------------------------------------
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)                  :: this
    INTEGER      ,INTENT   (IN),DIMENSION(:)        :: id
    INTEGER      ,INTENT   (IN)                     :: nl
    INTEGER      ,INTENT   (IN)                     :: ngvdb
    INTEGER                                         :: t,i,j,m
    REAL(KIND=GP),INTENT   (IN),DIMENSION(:)        :: lx,ly,lz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(3,ngvdb)  :: gvdb

    IF (this%myrank_.NE.0) THEN
!$omp parallel do
      DO i = 1,nl
        this%itop_(i) = i
      END DO 
!      CALL GPartComm_PPackV(this,this%stbuff_(:,1),this%nbuff_,id,lx,ly,lz,nl,this%itop_,nl)
      CALL GPartComm_PPackV2(this,this%stbuffp_,id,lx,ly,lz,this%itop_,nl)
      CALL GTStart(this%hcomm_)
!      CALL MPI_ISEND(this%stbuff_(:,1),4*nl+1,GC_REAL,0,this%myrank_,this%comm_,this%itsh_(1),this%ierr_)
      CALL MPI_ISEND(this%stbuffp_,nl,MPI_GPDataPackType,0,this%myrank_,this%comm_,this%itsh_(1),this%ierr_)
      CALL MPI_WAIT(this%itsh_(1),this%istatus_,this%ierr_)
      CALL GTAcc(this%hcomm_)
    ELSE
!$omp parallel do
      DO j = 1,nl
        i = id(j) + 1
        gvdb(1,i) = lx(j)
        gvdb(2,i) = ly(j)
        gvdb(3,i) = lz(j)
      END DO
    
      DO t = 1,this%nprocs_-1
        CALL GTStart(this%hcomm_)
        CALL MPI_RECV(this%rbbuffp_,this%maxparts_,MPI_GPDataPackType,t, &
                      t,this%comm_,this%ibrh_(1),this%ierr_)
        CALL MPI_GET_COUNT(this%ibrh_,MPI_GPDataPackType,m)
        CALL GTAcc(this%hcomm_)
!        m = int(this%rbbuff_(1,1))
!$omp parallel do
        DO j = 1,m
!          i = int(this%rbbuff_(4*j-2,1)) + 1
!          gvdb(1,i) = this%rbbuff_(4*j-1,1)
!          gvdb(2,i) = this%rbbuff_(4*j+0,1)
!          gvdb(3,i) = this%rbbuff_(4*j+1,1)
          i = this%rbbuffp_(j)%idp_ + 1
          gvdb(1,i) = this%rbbuffp_(j)%rp_(1)
          gvdb(2,i) = this%rbbuffp_(j)%rp_(2)
          gvdb(3,i) = this%rbbuffp_(j)%rp_(3)
        END DO
      END DO
    END IF

 END SUBROUTINE GPartComm_VDBSynch_t0
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_VDBSynch(this,gvdb,ngvdb,id,lx,ly,lz,nl,ptmp)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : VDBSynch
!  DESCRIPTION: Synch up global VDB from local vector data
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    gvdb    : global VDB containing part. position records, returned.
!    ngvdb   : no. records in global VDB. Fixed on entry.
!    id      : local part. ids
!    lx,ly,lz: local part. d.b. vectors
!    nl      : no. parts. in local pdb
!    ptmp    : tmp array of size of gvdb
!-----------------------------------------------------------------
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)                  :: this
    INTEGER      ,INTENT   (IN),DIMENSION(:)        :: id
    INTEGER      ,INTENT   (IN)                     :: nl
    INTEGER      ,INTENT   (IN)                     :: ngvdb
    INTEGER                                         :: i,j
    REAL(KIND=GP),INTENT   (IN),DIMENSION(:)        :: lx,ly,lz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(3,ngvdb)  :: gvdb,ptmp

!   CALL GTStart(this%hcomm_)
!   CALL MPI_ALLREDUCE(nl,ng,1,MPI_INTEGER,   &
!                      MPI_SUM,this%comm_,this%ierr_)
!   CALL GTAcc(this%hcomm_)

!!  IF ( this%myrank_.EQ.0 .AND. ng.NE.ngvdb ) THEN
!!    IF ( .NOT.present(scaller) ) THEN
!!      WRITE(*,*)'GPartComm_VDBSynch: inconsistent d.b.: expected: ', &
!!                 ngvdb, '; found: ',ng
!!    ELSE
!!      WRITE(*,*)'GPartComm_VDBSynch: caller:',trim(scaller),': inconsistent d.b.: expected: ', &
!!                 ngvdb, '; found: ',ng
!!    ENDIF
!!    STOP
!!  ENDIF

    DO j = 1, ngvdb
      gvdb(1:3,j) = 0.0_GP
      ptmp(1:3,j) = 0.0_GP
    ENDDO

    DO j = 1, nl
      i = id(j) + 1
      ptmp(1,i) = lx(j)
      ptmp(2,i) = ly(j)
      ptmp(3,i) = lz(j)
    ENDDO
    CALL GTStart(this%hcomm_)
    CALL MPI_ALLREDUCE(ptmp,gvdb,3*ngvdb,GC_REAL,   &
                       MPI_SUM,this%comm_,this%ierr_)
    CALL GTAcc(this%hcomm_)

 END SUBROUTINE GPartComm_VDBSynch
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_LagSynch_t0(this,gs,ngs,id,ls,nl)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : LagSynch
!  DESCRIPTION: Synch up task 0 Lagrangian scalar from local scalar data
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    gs      : global scalar containing 'synched' records, returned (task 0 only).
!    ngs     : no. records in global scalar. Fixed on entry.
!    id      : local part. ids
!    ls      : local scalar
!    nl      : no. parts. in local Lag. scalar
!-----------------------------------------------------------------
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)                  :: this
    INTEGER      ,INTENT   (IN),DIMENSION(:)        :: id
    INTEGER      ,INTENT   (IN)                     :: nl
    INTEGER      ,INTENT   (IN)                     :: ngs
    INTEGER                                         :: i,j,t,m
    REAL(KIND=GP),INTENT   (IN),DIMENSION(:)        :: ls
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:)        :: gs

    IF (this%myrank_.NE.0) THEN
      this%stbuff_(1,1) = REAL(nl,kind=GP)
!$omp parallel do
      DO i = 1,nl
        this%stbuffp_(i)%idp_   = id(i)
        this%stbuffp_(i)%rp_(1) = ls(i)
!        this%stbuff_(2*i,1) = REAL(id(i),kind=GP)
!        this%stbuff_(2*i+1,1) = ls(i)
      END DO 
      CALL GTStart(this%hcomm_)
!      CALL MPI_ISEND(this%stbuff_,2*nl+1,GC_REAL,0,this%myrank_,this%comm_,this%ibsh_(1),this%ierr_)
      CALL MPI_ISEND(this%stbuffp_,nl,MPI_GPDataPackType,0,this%myrank_,this%comm_,this%itsh_(1),this%ierr_)
      CALL MPI_WAIT(this%itsh_(1),this%istatus_,this%ierr_)
      CALL GTAcc(this%hcomm_)
    ELSE
!$omp parallel do
      DO j = 1,nl
        i = id(j) + 1
        gs(i) = ls(j)
      END DO
     
      DO t = 1,this%nprocs_-1
        CALL GTStart(this%hcomm_)
        CALL MPI_RECV(this%rbbuffp_,this%maxparts_,MPI_GPDataPackType,t, &
                      t,this%comm_,this%ibrh_(1),this%ierr_)
        CALL MPI_GET_COUNT(this%ibrh_,MPI_GPDataPackType,m)
!        CALL MPI_WAIT(this%ibrh_(1),this%istatus_,this%ierr_)
        CALL GTAcc(this%hcomm_)
!        m = int(this%rbbuff_(1,1))
!$omp parallel do
        DO j = 1,m
!          i = int(this%rbbuff_(2*j,1)) + 1
!          gs(i) = this%rbbuff_(2*j+1,1)
          i     = this%rbbuffp_(j)%idp_ + 1
          gs(i) = this%rbbuffp_(j)%rp_(1)
        END DO
      END DO
    END IF

  END SUBROUTINE GPartComm_LagSynch_t0
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPartComm_LagSynch(this,gs,ngs,id,ls,nl,ptmp)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : LagSynch
!  DESCRIPTION: Synch up global Lagrangian scalar from local scalar data
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    gs      : global scalar containing 'synched' records, returned.
!    ngs     : no. records in global scalar. Fixed on entry.
!    id      : local part. ids
!    ls      : local scalar
!    nl      : no. parts. in local Lag. scalar
!    ptmp    : tmp array of size of gs
!-----------------------------------------------------------------
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)                  :: this
    INTEGER      ,INTENT   (IN),DIMENSION(*)        :: id
    INTEGER      ,INTENT   (IN)                     :: nl
    INTEGER      ,INTENT   (IN)                     :: ngs
    INTEGER                                         :: i,j
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)        :: ls
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)        :: gs,ptmp


    DO j = 1, ngs
      gs  (j) = 0.0_GP
      ptmp(j) = 0.0_GP
    ENDDO

    DO j = 1, nl
      i = id(j) + 1
      ptmp(i) = ls(j)
    ENDDO
    CALL GTStart(this%hcomm_)
    CALL MPI_ALLREDUCE(ptmp,gs,ngs,GC_REAL,   &
                       MPI_SUM,this%comm_,this%ierr_)
    CALL GTAcc(this%hcomm_)

  END SUBROUTINE GPartComm_LagSynch
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_SetCacheParam(this,csize,nstrip)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : SetCacheParam
!  DESCRIPTION: Set cache size and strip-mining size for transpose
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    csize   : cache-size
!    nstrip  : strip mining size
!-----------------------------------------------------------------
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)                  :: this
    INTEGER      ,INTENT   (IN)                     :: csize,nstrip

    this%csize_  = csize
    this%nstrip_ = nstrip

  END SUBROUTINE GPartComm_SetCacheParam
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_Copy2Ext(this,vext,v)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Copy2Ext
!  DESCRIPTION: Copy field from regular to extended grid
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    vext    : extended-grid field
!    v       : regular-grid field
!    ldims   : local dims of v
!-----------------------------------------------------------------
!$  USE threads
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)                      :: this
    INTEGER                                             :: i,j,jm,k,km,ngp,ngz,nex,nexy
    INTEGER                                             :: nx,nxy,ny
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: v
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: vext

    ngz  = this%nzghost_
    ngp  = ngz * this%iextperp_
    nexy = (this%nd_(1)+2*ngp) * (this%nd_(2)+2*ngp)
    nex  = this%nd_(1)+2*ngp
    nx   = this%nd_(1)
    ny   = this%nd_(2)
    nxy  = nx*ny

!$omp parallel do if(ny*(this%kend_-this%ksta_+1).GE.nth) collapse(2) private(jm,i)
     DO km = 0,this%kend_-this%ksta_
!    DO k = 1,this%kend_-this%ksta_+1
!      km = k-1
      DO j=1,ny
        jm = j-1
        DO i=1,nx
          vext(i+ngp+(jm+ngp)*nex+(km+ngz)*nexy) = v(i+jm*nx+km*nxy)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE GPartComm_Copy2Ext
!-----------------------------------------------------------------
!-----------------------------------------------------------------

 SUBROUTINE GPartComm_ResizeArrays(this,newmparts,onlyinc)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Resize_Arrays
!  DESCRIPTION: Resize all arrays in the GPart class (including 
!               subclases, i.e. communicator, spline)
!  ARGUMENTS  :
!    this     : 'this' class instance
!    newmparts: new number of particles
!    onlyinc  : if true, will only resize to increase array size
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)  :: this
    INTEGER         ,INTENT(IN)     :: newmparts
    LOGICAL         ,INTENT(IN)     :: onlyinc
    INTEGER                         :: n,nzg,nd(3),newbuff

    nzg = this%nzghost_
    nd  = this%nd_

    this%maxparts_ = newmparts

    n = SIZE(this%ibot_)
    IF ((n.lt.newmparts).OR.((n.gt.newmparts).AND..NOT.onlyinc)) THEN
      CALL Resize_IntArray(this%ibot_ ,newmparts,.true. )
      CALL Resize_IntArray(this%itop_ ,newmparts,.true.)
      CALL Resize_IntArray(this%oldid_,newmparts,.false.)
      CALL Resize_DataPack(this%sbbuffp_,newmparts)
      CALL Resize_DataPack(this%stbuffp_,newmparts)
      CALL Resize_DataPack(this%rbbuffp_,newmparts)
      CALL Resize_DataPack(this%rtbuffp_,newmparts)
      this%maxparts_ = newmparts
      IF ( this%intrfc_ .GE. 1 ) THEN
!        newbuff = MAX(newmparts*4,3*nd(1)*nd(2)*nzg+nzg+1) + 1
        newbuff = 3*nd(1)*nd(2)*nzg+nzg+2
      ELSE
!        newbuff = MAX(newmparts*4,  nd(1)*nd(2)*nzg+nzg+1) + 1
        newbuff =   nd(1)*nd(2)*nzg+nzg+2
      ENDIF
      n = this%nbuff_
      IF ((n.lt.newbuff).OR.((n.gt.newbuff).AND..NOT.onlyinc)) THEN
        CALL Resize_ArrayRank2Transposed(this%sbbuff_,newbuff,.false.)
        CALL Resize_ArrayRank2Transposed(this%stbuff_,newbuff,.false.)
        CALL Resize_ArrayRank2Transposed(this%rbbuff_,newbuff,.false.)
        CALL Resize_ArrayRank2Transposed(this%rtbuff_,newbuff,.false.)
        this%nbuff_ = newbuff
      ENDIF
    ENDIF

    RETURN
  END SUBROUTINE GPartComm_ResizeArrays
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE Resize_DataPack(a,new_size)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Resize_DataPack
!  DESCRIPTION: Resize input 1D GPData array to new size, without copying
!  ARGUMENTS  :
!    a       : array to be resized
!    new_size: new size for array
!-----------------------------------------------------------------
!$  USE threads

    IMPLICIT NONE
    TYPE(GPData),INTENT(INOUT),ALLOCATABLE,DIMENSION(:) :: a
    INTEGER,INTENT(IN)                                  :: new_size
    TYPE(GPData)              ,ALLOCATABLE,DIMENSION(:) :: temp

    ALLOCATE ( temp(new_size) )

    CALL move_alloc(temp, a)

    RETURN
  END SUBROUTINE Resize_DataPack
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE Resize_IntArray(a,new_size,docpy)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Resize_ArrayRank1
!  DESCRIPTION: Resize input 1D integer array to new size
!  ARGUMENTS  :
!    a       : array to be resized
!    new_size: new size for array
!    docpy   : if true, keeps previous data in resized array
!-----------------------------------------------------------------
!$  USE threads

    IMPLICIT NONE
    INTEGER,INTENT(INOUT),ALLOCATABLE,DIMENSION(:) :: a
    INTEGER,INTENT(IN)                             :: new_size
    LOGICAL,INTENT(IN)                             :: docpy
    INTEGER              ,ALLOCATABLE,DIMENSION(:) :: temp
    INTEGER                                        :: i,n

    ALLOCATE ( temp(new_size) )
    IF (docpy) THEN
      n = SIZE(a)
      n = MIN(n,new_size)
!$omp parallel do if(n.gt.NMIN_OMP)
      DO i = 1,n
        temp(i) = a(i)
      END DO
    END IF

    CALL move_alloc(temp, a)

    RETURN
  END SUBROUTINE Resize_IntArray
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE Resize_IntArrayRank2(a,new_size,docpy)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Resize_IntArrayRank2
!  DESCRIPTION: Resize input integer 2D array to new size
!  ARGUMENTS  :
!    a       : array to be resized
!    new_size: new size for array (second axis)
!    docpy   : if true, keeps previous data in resized array
!-----------------------------------------------------------------
!$  USE threads

    IMPLICIT NONE
    INTEGER,INTENT(INOUT),ALLOCATABLE,DIMENSION(:,:) :: a
    INTEGER      ,INTENT(IN)                         :: new_size
    LOGICAL      ,INTENT(IN)                         :: docpy
    INTEGER              ,ALLOCATABLE,DIMENSION(:,:) :: temp
    INTEGER                                          :: i,n,shp(2)

    shp = SHAPE(a)
    ALLOCATE ( temp(shp(1),new_size) )

    IF (docpy) THEN
      n = shp(2)
      n = MIN(n,new_size)
!$omp parallel do if(n.gt.NMIN_OMP)
      DO i = 1,n
        temp(:,i) = a(:,i)
      END DO
    END IF

    CALL move_alloc(temp, a)

    RETURN
  END SUBROUTINE Resize_IntArrayRank2
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE Resize_ArrayRank1(a,new_size,docpy)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Resize_ArrayRank1
!  DESCRIPTION: Resize input 1D array to new size
!  ARGUMENTS  :
!    a       : array to be resized
!    new_size: new size for array
!    docpy   : if true, keeps previous data in resized array
!-----------------------------------------------------------------
!$  USE threads

    IMPLICIT NONE
    REAL(KIND=GP),INTENT(INOUT),ALLOCATABLE,DIMENSION(:) :: a
    INTEGER      ,INTENT(IN)                             :: new_size
    LOGICAL      ,INTENT(IN)                             :: docpy
    REAL(KIND=GP)              ,ALLOCATABLE,DIMENSION(:) :: temp
    INTEGER                                              :: i,n

    ALLOCATE ( temp(new_size) )
    IF (docpy) THEN
      n = SIZE(a)
      n = MIN(n,new_size)
!$omp parallel do if(n.gt.NMIN_OMP)
      DO i = 1,n
        temp(i) = a(i)
      END DO
    END IF

    CALL move_alloc(temp, a)

    RETURN
  END SUBROUTINE Resize_ArrayRank1
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE Resize_ArrayRank2(a,new_size,docpy)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Resize_ArrayRank2
!  DESCRIPTION: Resize input 2D array to new size
!  ARGUMENTS  :
!    a       : array to be resized
!    new_size: new size for array (second axis)
!    docpy   : if true, keeps previous data in resized array
!-----------------------------------------------------------------
!$  USE threads

    IMPLICIT NONE
    REAL(KIND=GP),INTENT(INOUT),ALLOCATABLE,DIMENSION(:,:) :: a
    INTEGER      ,INTENT(IN)                               :: new_size
    LOGICAL      ,INTENT(IN)                               :: docpy
    REAL(KIND=GP)              ,ALLOCATABLE,DIMENSION(:,:) :: temp
    INTEGER                                                :: i,n,shp(2)

    shp = SHAPE(a)
    ALLOCATE ( temp(shp(1),new_size) )

    IF (docpy) THEN
      n = shp(2)
      n = MIN(n,new_size)
!$omp parallel do if(n.gt.NMIN_OMP)
      DO i = 1,n
        temp(:,i) = a(:,i)
      END DO
    END IF

    CALL move_alloc(temp, a)

    RETURN
  END SUBROUTINE Resize_ArrayRank2
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE Resize_ArrayRank2Transposed(a,new_size,docpy)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Resize_ArrayRank2Transposed
!  DESCRIPTION: Resize input 2D array to new size
!  ARGUMENTS  :
!    a       : array to be resized
!    new_size: new size for array (first axis)
!    docpy   : if true, keeps previous data in resized array
!-----------------------------------------------------------------
!$  USE threads
 
    IMPLICIT NONE
    REAL(KIND=GP),INTENT(INOUT),ALLOCATABLE,DIMENSION(:,:) :: a
    INTEGER      ,INTENT(IN)                               :: new_size
    LOGICAL      ,INTENT(IN)                               :: docpy
    REAL(KIND=GP)              ,ALLOCATABLE,DIMENSION(:,:) :: temp
    INTEGER                                                :: i,n,j,shp(2)

    shp = SHAPE(a) 
    ALLOCATE ( temp(new_size,shp(2)) )
    
    IF (docpy) THEN
      n = shp(1)
      n = MIN(n,new_size)
!$omp parallel do if((n*shp(2)).gt.NMIN_OMP) collapse(2)
      DO i = 1,n
        DO j = 1,shp(2)
          temp(i,j) = a(i,j)
        END DO
      END DO
    END IF
    
    CALL move_alloc(temp, a)

    RETURN
  END SUBROUTINE Resize_ArrayRank2Transposed
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!  INCLUDE 'gfieldcomm_contain.f90'
!!! GFIELDCOMM SUBROUTINES UNTIL THE END
  SUBROUTINE GPartComm_SlabDataReturnSF(this,v,vext,method)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPartComm_SlabDataReturnSF
!  DESCRIPTION: Does bdy exchange of field component, v. Output
!               is to data on regular grids, v. 'SF' means
!               that this is the 'single-field' interface.
!  ARGUMENTS  :
!    this      : 'this' class instance (IN)
!    v         : Eulerian velocity component returned on regular
!                grid.
!    vext      : Eulerian velocity components on extended grid.
!    method    : Ghost zone aggregation method
!                   = UNPACK_REP for replacement
!                   = UNPACK_SUM for sum
!
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)                      :: this
    INTEGER      ,INTENT   (IN)                         :: method
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: v
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: vext

    INTEGER                                             :: itask,j,buffsize

    IF ( this%intrfc_ .GE. 1 ) THEN
      WRITE(*,*) 'GPartComm_SlabDataReturnSF: SF interface expected'
      STOP
    ENDIF
    IF ( this%nprocs_ .EQ. 1 ) THEN
      CALL GPartComm_LocalDataRetSF(this,vext,v,method)
      RETURN
    ENDIF

    buffsize = this%nd_(1)*this%nd_(2)*this%nzghost_ + this%nzghost_ + 2

    CALL GPartComm_Copy2Reg(this,v,vext)

    ! post receives:
    CALL GTStart(this%hcomm_)
    DO j=1,this%nbsnd_  ! from bottom task:
      itask = this%ibsndp_(j)
      CALL MPI_IRECV(this%rbbuff_(:,j),2*this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%ibrh_(j),this%ierr_)
    ENDDO
    CALL GTAcc(this%hcomm_)

    ! return data:
    DO j=1,this%nbret_  ! to bottom task:
      itask = this%ibretp_(j)
!      PRINT *, 'Sending', j, this%myrank_, itask
      CALL GPartComm_PackRetSF(this,this%sbbuff_(:,j),vext,j,'b')
      CALL GTStart(this%hcomm_)
      CALL MPI_ISEND(this%sbbuff_,buffsize,GC_REAL,itask, &
                     1,this%comm_,this%ibsh_(j),this%ierr_)
      CALL GTAcc(this%hcomm_)
    ENDDO
!

    CALL GTStart(this%hcomm_)
    DO j=1,this%ntsnd_  ! from top task:
      itask = this%itsndp_(j)
      CALL MPI_IRECV(this%rtbuff_(:,j),2*this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%itrh_(j),this%ierr_)
    ENDDO
    CALL GTAcc(this%hcomm_)


    DO j=1,this%ntret_  ! to top task:
      itask = this%itretp_(j)
      CALL GPartComm_PackRetSF(this,this%stbuff_(:,j),vext,j,'t')
      CALL GTStart(this%hcomm_)
      CALL MPI_ISEND(this%stbuff_,buffsize,GC_REAL,itask, &
                     1,this%comm_,this%itsh_(j),this%ierr_)
      CALL GTAcc(this%hcomm_)
    ENDDO

    CALL GTStart(this%hcomm_)
    DO j=1,this%nbret_
      CALL MPI_WAIT(this%ibsh_(j),this%istatus_,this%ierr_)
    ENDDO
    DO j=1,this%nbsnd_
      CALL MPI_WAIT(this%ibrh_(j),this%istatus_,this%ierr_)
    ENDDO


    DO j=1,this%ntret_
      CALL MPI_WAIT(this%itsh_(j),this%istatus_,this%ierr_)
    ENDDO
    DO j=1,this%ntsnd_
      CALL MPI_WAIT(this%itrh_(j),this%istatus_,this%ierr_)
    ENDDO
    CALL GTAcc(this%hcomm_)


    ! Unpack received data:
    DO j=1,this%nbsnd_
       CALL GPartComm_UnpackRetSF(this,v,this%rbbuff_(:,j),&
                          this%nbuff_,'b',method,this%ierr_)
    ENDDO
    DO j=1,this%ntsnd_
       CALL GPartComm_UnpackRetSF(this,v,this%rtbuff_(:,j),&
                          this%nbuff_,'t',method,this%ierr_)
    ENDDO

  END SUBROUTINE GPartComm_SlabDataReturnSF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_PackRetSF(this,buff,vext,isnd,sdir)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PPackSF
!  DESCRIPTION: packs ret buffer with (single) field
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    buff    : packed buffer (returned)
!    vext    : Eulerian velocity component on extended grid
!              in phys. space (IN)
!    isnd    : which send this is
!    sdir    : 't' for top, 'b' for bottom
!
!-----------------------------------------------------------------
!$  USE threads
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)          :: this
    INTEGER      ,INTENT   (IN)             :: isnd
    INTEGER                                 :: i,j,k,m,nt,nx,ny,nxy
    INTEGER                                 :: jm,km
    REAL(KIND=GP),INTENT  (OUT),DIMENSION(*):: buff
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*):: vext
    CHARACTER*(*),INTENT   (IN)             :: sdir


    IF ( sdir(1:1).NE.'b' .AND. sdir(1:1).NE.'B' &
    .AND.sdir(1:1).NE.'t' .AND. sdir(1:1).NE.'T' ) THEN
      WRITE(*,*) 'GPartComm_PackRetSF: Bad direction descriptor'
      STOP
    ENDIF

    nx  = this%nd_(1)
    ny  = this%nd_(2)
    nxy = nx*ny
    IF      ( sdir(1:1) .EQ. 'b' .OR. sdir(1:1) .EQ. 'B' ) THEN
    ! Pack for send to rank at bottom:
    !  ...header
      nt = 1
      buff(1)  = this%ibretnz_(isnd)      ! no. z-indices included
      DO j = 1, this%ibretnz_(isnd)
        nt       = nt + 1
        buff(nt) = this%ibretdst_(isnd,j) ! z-index in regular grid
      ENDDO

    !  ...data
      DO m = 1,this%ibretnz_(isnd)
        k = this%ibret_(isnd,m)
        km = k-1
!$omp parallel do if(ny.ge.nth) private(jm,i)
        DO j = 1,ny
          jm = j-1
          DO i = 1,nx
!            nt = nt + 1
            buff(nt+jm*nx+i) = vext(i+jm*nx+km*nxy)
          ENDDO
        ENDDO
        nt = nt + nxy
      ENDDO

    ELSE !  Pack for send to rank at top:

      ! ...header
      nt = 1
      buff(1)  = this%itretnz_(isnd)      ! no. z-indices included
      DO j = 1, this%itretnz_(isnd)
        nt       = nt + 1
        buff(nt) = this%itretdst_(isnd,j) ! z-index in regular grid
      ENDDO

      ! ...data
      DO m = 1,this%itretnz_(isnd)
        k = this%itret_(isnd,m)
        km = k-1
!$omp parallel do if(ny.ge.nth) private(jm,i)
        DO j = 1,ny
          jm = j-1
          DO i = 1,nx
!            nt = nt + 1
            buff(nt+jm*nx+i) = vext(i+jm*nx+km*nxy)
          ENDDO
        ENDDO
        nt = nt + nxy
      ENDDO

    ENDIF

  END SUBROUTINE GPartComm_PackRetSF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_UnpackRetSF(this,v,buff,nbuff,sb,method,ierr)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : RetUnpackSF
!  DESCRIPTION: Unpacks recv buffer with into regular (single) field.
!               Messages received are self-referential, so contain info
!               on where to 'return' recvd data. So, there is no 't'
!               or 'b' designation required for unpacking.
!  ARGUMENTS  :
!    this        : 'this' class instance (IN)
!    v           : Eulerian velocity component on regular grid
!                  in phys. space (IN)
!    buff        : packed buffer (input) from which to store into
!                  extended grid quantities.
!    nbuff       : maximum buff size
!    sb          : optional buffer name ,'t' or 'b'.
!    method      : data aggregation method
!                   = UNPACK_REP for replacement
!                   = UNPACK_SUM for sum
!    ierr        : err flag: 0 if success; else 1
!
!-----------------------------------------------------------------
!$  USE threads
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)          :: this
    INTEGER                                 :: i,j,k,m,ngp,nx,nex,nxy,nexy,ny,nz
    INTEGER                                 :: im,ip,ir,ixy,iz,jm,km,nh
    INTEGER      ,INTENT   (IN)             :: nbuff,method
    INTEGER      ,INTENT(INOUT)             :: ierr ! not used now
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*):: buff
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*):: v
    CHARACTER(len=1),INTENT(IN),OPTIONAL    :: sb

    ierr = 0
    nx   = this%nd_(1)
    ny   = this%nd_(2)
    nxy  = nx*ny
    ngp  = this%nzghost_*this%iextperp_
    nex  = nx+2*ngp
    nexy = (nx+2*ngp)*(ny+2*ngp)

    ! Unpack from either buffer:
    ! For each task, message is of form:
    !     #z-indices:z-index_0:z-index_1:...:nx*ny_0:nx*ny_1: ...
    nz = int(buff(1))
    nh = nz + 1 ! no. items in header
    IF ( method .EQ. UNPACK_REP ) THEN
      DO m = 1,nz
        k   = int(buff(m+1))
        km  = k-1
!        ixy = 1
!$omp parallel do if(ny.ge.nth) private(jm,i,ixy,im,ir)
        DO j = 1,ny
          jm = j-1
          DO i = 1,nx
            ixy = jm*nx + i
            im = i+ngp+(jm+ngp)*nex+km*nexy
            ir = nh+(m-1)*nxy+ixy
            v(im) = buff(ir)
!            ixy = ixy + 1
          ENDDO
        ENDDO
      ENDDO
    ELSE IF ( method .EQ. UNPACK_SUM ) THEN
      DO m = 1,nz
        k   = int(buff(m+1))
        km  = k-1
!        ixy = 1
!$omp parallel do if(ny.ge.nth) private(jm,i,ixy,im,ir)
        DO j = 1,ny
          jm = j-1
          DO i = 1,nx
            ixy = jm*nx + i
            im = i+ngp+(jm+ngp)*nex+km*nexy
            ir = nh+(m-1)*nxy+ixy
            v(im) = v(im) + buff(ir)
!            ixy = ixy + 1
          ENDDO
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE GPartComm_UnpackRetSF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_LocalDataRetSF(this,vext,v,method)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : LocalDataExchSF
!  DESCRIPTION: Does 'bdy exchange' of (single) velocity component, when there's
!               only a single MPI task. This is a single-field interface.
!  ARGUMENTS  :
!    this              : 'this' class instance (IN)
!    v                 : Eulerian velocity components returned on regular grid.
!                        Must be of size nd_ set in constructor.
!    vext              : Eulerian velocity components on extended grid (that
!                        used to hold ghost zones).
!    method            : data aggregation method
!                         = UNPACK_REP for replacement
!                         = UNPACK_SUM for sum
!
!-----------------------------------------------------------------
!$  USE threads
    USE mpivars

    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)                      :: this
    INTEGER                                             :: i,j,k,ngp,ngz,nex,nexy,nez
    INTEGER                                             :: nx,nxy,ny,nz
    INTEGER      ,INTENT   (IN)                         :: method
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: v
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: vext

    ngz  = this%nzghost_
    ngp  = ngz * this%iextperp_
    nexy = (this%nd_(1)+2*ngp) * (this%nd_(2)+2*ngp)
    nex  = this%nd_(1)+2*ngp
    nez  = this%nd_(3)+  ngz
    nx   = this%nd_(1)
    ny   = this%nd_(2)
    nz   = this%nd_(3)
    nxy  = nx*ny

    CALL GPartComm_Copy2Reg(this,v,vext)
    IF ( method .EQ. UNPACK_REP ) THEN
!$omp parallel do if(ngz*ny.ge.nth) collapse(2) private(i)
      DO k = 1,ngz  ! extended zones
        DO j = 1,ny
          DO i = 1,nx
            ! set top bcs:
            v(i+(j-1)*nx+       (k-1)*nxy) = vext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy)

            ! set bottom bcs:
            v(i+(j-1)*nx+(nz-ngz+k-1)*nxy) = vext(i+ngp+(j+ngp-1)*nex+    (k-1)*nexy)

          ENDDO
        ENDDO
      ENDDO
    ELSE IF ( method .EQ. UNPACK_SUM ) THEN
!$omp parallel do if(ngz*ny.ge.nth) collapse(2) private(i)
      DO k = 1,ngz  ! extended zones
        DO j = 1,ny
          DO i = 1,nx
            ! set top bcs:
            v(i+(j-1)*nx+(k-1)*nxy)        = &
                   v(i+(j-1)*nx+(k-1)*nxy) + vext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy)

            ! set bottom bcs:
            v(i+(j-1)*nx+(nz-ngz+k-1)*nxy) = &
                v(i+(j-1)*nx+(nz-ngz+k-1)*nxy) + vext(i+ngp+(j+ngp-1)*nex+(k-1)*nexy)

          ENDDO
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE GPartComm_LocalDataRetSF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_SlabDataReturnMF(this,vx,vy,vz,vxext,vyext,vzext,method)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPartComm_SlabDataReturnMF
!  DESCRIPTION: Does bdy exchange of velocity component, vx,vy,vz. Output
!               is to data on extended grids, vxext, vyext, vzexy. 'MF'
!               means that this is the 'multi-field' interface.
!  ARGUMENTS  :
!    this              : 'this' class instance (IN)
!    vx,vy,vz          : Eulerian velocity components returned on regular
!                        grid
!    vxext,vyext,vzext : Eulerian velocity components on extended grid.
!    method            : Ghost zone aggregation method
!                         = UNPACK_REP for replacement
!                         = UNPACK_SUM for sum
!
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)                      :: this
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: vxext,vyext,vzext
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: vx,vy,vz
    INTEGER      ,INTENT   (IN)                         :: method
    INTEGER                                             :: itask,j

    IF ( this%intrfc_ .LT. 1 ) THEN
      WRITE(*,*) 'GPartComm_SlabDataReturnMF: SF interface expected'
      STOP
    ENDIF

    IF ( this%nprocs_ .EQ. 1 ) THEN
      CALL GPartComm_LocalDataRetMF(this,vxext,vyext,vzext,vx,vy,vz,method)
      RETURN
    ENDIF

    ! Post receives:
    CALL GTStart(this%hcomm_)
    DO j=1,this%nbsnd_  ! from bottom task:
      itask = this%ibsndp_(j)
      CALL MPI_IRECV(this%rbbuff_(:,j),this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%ibrh_(j),this%ierr_)
    ENDDO
    CALL GTAcc(this%hcomm_)

    CALL GTStart(this%hcomm_)
    DO j=1,this%ntsnd_  ! from top task:
      itask = this%itsndp_(j)
      CALL MPI_IRECV(this%rtbuff_(:,j),this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%itrh_(j),this%ierr_)
    ENDDO
    CALL GTAcc(this%hcomm_)

    !
    ! return data:
    DO j=1,this%nbret_  ! to bottom task:
      itask = this%ibretp_(j)
      CALL GPartComm_PackRetMF(this,this%sbbuff_(:,j),vxext,vyext,vzext,j,'b')
      CALL GTStart(this%hcomm_)
      CALL MPI_ISEND(this%sbbuff_,this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%ibsh_(j),this%ierr_)
      CALL GTAcc(this%hcomm_)
    ENDDO
    DO j=1,this%ntret_  ! to top task:
      itask = this%itretp_(j)
      CALL GPartComm_PackRetMF(this,this%stbuff_(:,j),vxext,vyext,vzext,j,'t')
      CALL GTStart(this%hcomm_)
      CALL MPI_ISEND(this%stbuff_,this%nbuff_,GC_REAL,itask, &
                     1,this%comm_,this%itsh_(j),this%ierr_)
      CALL GTAcc(this%hcomm_)

    ENDDO

    CALL GTStart(this%hcomm_)
    DO j=1,this%nbret_
      CALL MPI_WAIT(this%ibsh_(j),this%istatus_,this%ierr_)
    ENDDO
    DO j=1,this%ntret_
      CALL MPI_WAIT(this%itsh_(j),this%istatus_,this%ierr_)
    ENDDO
    DO j=1,this%nbsnd_
      CALL MPI_WAIT(this%ibrh_(j),this%istatus_,this%ierr_)
    ENDDO
    DO j=1,this%ntsnd_
      CALL MPI_WAIT(this%itrh_(j),this%istatus_,this%ierr_)
    ENDDO
    CALL GTAcc(this%hcomm_)

    ! Unpack received data:
    DO j=1,this%nbsnd_
      CALL GPartComm_UnpackRetMF(this,vx,vy,vz,this%rbbuff_(:,j),method)
    ENDDO
    DO j=1,this%ntsnd_
      CALL GPartComm_UnpackRetMF(this,vx,vy,vz,this%rtbuff_(:,j),method)
    ENDDO

  END SUBROUTINE GPartComm_SlabDataReturnMF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_PackRetMF(this,buff,vxext,vyext,vzext,isnd,sdir)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PackMF
!  DESCRIPTION: packs ret buffer with fields; multi-field interface
!  ARGUMENTS  :
!    this             : 'this' class instance (IN)
!    buff             : packed buffer (returned)
!    vxext,vyext,vzext: Eulerian velocity component on extended
!                       grid in phys. space (IN)
!    isnd             : which send this is
!    sdir             : 't' for top, 'b' for bottom
!
!-----------------------------------------------------------------
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)          :: this
    INTEGER      ,INTENT   (IN)             :: isnd
    INTEGER                                 :: i,j,k,m,nt,nx,nxy,ny
    INTEGER                                 :: jm,km
    REAL(KIND=GP),INTENT  (OUT),DIMENSION(*):: buff
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*):: vxext,vyext,vzext
    CHARACTER*(*),INTENT   (IN)             :: sdir


    IF ( sdir(1:1).NE.'b' .AND. sdir(1:1).NE.'B' &
    .AND.sdir(1:1).NE.'t' .AND. sdir(1:1).NE.'T' ) THEN
      WRITE(*,*) 'GPartComm_PackRetMF: Bad direction descriptor'
      STOP
    ENDIF

    nx  = this%nd_(1)
    ny  = this%nd_(2)
    nxy = nx*ny
    IF      ( sdir(1:1) .EQ. 'b' .OR. sdir(1:1) .EQ. 'B' ) THEN
    ! Pack for send to rank at bottom:
    !  ...header
      nt = 1
      buff(1)  = this%ibretnz_(isnd)       ! no. z-indices included
      DO j = 1, this%ibretnz_(isnd)
        nt       = nt + 1
        buff(nt) = this%ibretdst_(isnd,j) ! z-index in regular grid
      ENDDO

    !  ...data
      DO m = 1,this%ibretnz_(isnd)
        k = this%ibret_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            buff(nt) = vxext(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

      DO m = 1, this%ibretnz_(isnd)
        k = this%ibret_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            buff(nt) = vyext(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

      DO m = 1,this%ibretnz_(isnd)
        k = this%ibret_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            buff(nt) = vzext(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

    ELSE !  Pack for send to rank at top:

      ! ...header
      nt = 1
      buff(1)  = this%itretnz_(isnd)      ! no. z-indices included
      DO j = 1, this%itretnz_(isnd)
        nt       = nt + 1
        buff(nt) = this%itretdst_(isnd,j) ! z-index in regular grid
      ENDDO

      ! ...data
      DO m = 1,this%itretnz_(isnd)
        k = this%itret_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            buff(nt) = vxext(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

      DO m = 1,this%itretnz_(isnd)
        k = this%itret_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            buff(nt) = vyext(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

      DO m = 1,this%itretnz_(isnd)
        k = this%itret_(isnd,m)
        km = k-1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            buff(nt) = vzext(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE GPartComm_PackRetMF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_UnpackRetMF(this,vx,vy,vz,buff,method)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : UnpackMF
!  DESCRIPTION: Unpacks recv buffer with into extended (single) field
!               Messages received are self-referential, so contain info
!               on where to 'send' recvd data. So, there is no 't' or
!               'b' designation required for unpacking.
!  ARGUMENTS  :
!    this        : 'this' class instance (IN)
!    buff        : packed buffer (input) from which to store into
!                  extended grid quantities.
!    vx,vy,vz    : Eulerian velocity component on regular grid
!                  in phys. space (IN)
!    method      : data aggregation method
!                   = UNPACK_REP for replacement
!                   = UNPACK_SUM for sum
!
!-----------------------------------------------------------------
    USE mpivars
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)          :: this
    INTEGER                                 :: i,j,k,m,ngp,ngz,nx,nxy,ny,nz
    INTEGER                                 :: ixy,jm,km,nh
    INTEGER      ,INTENT   (IN)             :: method
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*):: buff
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*):: vx,vy,vz

    nx  = this%nd_(1)
    ny  = this%nd_(2)
    nxy = nx*ny
    ngz = this%nzghost_;
    ngp = ngz*this%iextperp_

  ! Unpack from either top or bottom buffer:
    nz = int(buff(1))
    nh = nz + 1 ! no. items in header
    IF ( method .EQ. UNPACK_REP ) THEN
      DO m = 1,nz
        k   = int(buff(m+1))
        km  = k-1

        ixy = 1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            vx(i+ngp+(jm+ngp)*nx+km*nxy) = buff(nh+(m-1)*nxy+ixy)
            ixy = ixy + 1
          ENDDO
        ENDDO

        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            vy(i+ngp+(jm+ngp)*nx+km*nxy) = buff(nh+(m-1)*nxy+ixy)
            ixy = ixy + 1
          ENDDO
        ENDDO

        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            vz(i+ngp+(jm+ngp)*nx+km*nxy) = buff(nh+(m-1)*nxy+ixy)
            ixy = ixy + 1
          ENDDO
        ENDDO

      ENDDO
    ELSE IF ( method .EQ. UNPACK_SUM ) THEN
      DO m = 1,nz
        k   = int(buff(m+1))
        km  = k-1

        ixy = 1
        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            vx(i+ngp+(jm+ngp)*nx+km*nxy) = vx(i+ngp+(jm+ngp)*nx+km*nxy) &
                                          + buff(nh+(m-1)*nxy+ixy)
            ixy = ixy + 1
          ENDDO
        ENDDO

        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            vy(i+ngp+(jm+ngp)*nx+km*nxy) = vy(i+ngp+(jm+ngp)*nx+km*nxy) &
                                          + buff(nh+(m-1)*nxy+ixy)
            ixy = ixy + 1
          ENDDO
        ENDDO

        DO j = 1, ny
          jm = j-1
          DO i = 1, nx
            vz(i+ngp+(jm+ngp)*nx+km*nxy) = vz(i+ngp+(jm+ngp)*nx+km*nxy) &
                                          + buff(nh+(m-1)*nxy+ixy)
            ixy = ixy + 1
          ENDDO
        ENDDO

      ENDDO
    ENDIF

  END SUBROUTINE GPartComm_UnpackRetMF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_LocalDataRetMF(this,vxext,vyext,vzext,vx,vy,vz,method)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : LocalDataExch
!  DESCRIPTION: Does 'bdy exchange' of velocity component, when there's
!               only a single MPI task.
!  ARGUMENTS  :
!    this              : 'this' class instance (IN)
!    vxext,vyext,vzext : Eulerian velocity components returned on extended
!                        grid (that used to hold ghost zones). Only z-conditions
!                        are imposed; lateral periodicity is not handled here.
!                        Lateral ghost zones can be accounted for by setting
!                        this%iextperp_=1 in contructor.
!    vx,vy,vz          : Eulerian velocity components returned on regular
!                        grid. Must be of size nd_ set in constructor
!    method            : data aggregation method
!                         = UNPACK_REP for replacement
!                         = UNPACK_SUM for sum
!
!-----------------------------------------------------------------
    USE mpivars

    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)                      :: this
    INTEGER                                             :: i,j,k,ngp,ngz,nex,nexy,nez
    INTEGER                                             :: nx,nxy,ny,nz
    INTEGER                                             :: jm,km
    INTEGER      , INTENT  (IN)                         :: method
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: vx,vy,vz
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: vxext,vyext,vzext

    ngz  = this%nzghost_
    ngp  = ngz * this%iextperp_
    nexy = (this%nd_(1)+2*ngp) * (this%nd_(2)+2*ngp)
    nex  = this%nd_(1)+2*ngp
    nez  = this%nd_(3)+  ngz
    nx   = this%nd_(1)
    ny   = this%nd_(2)
    nz   = this%nd_(3)
    nxy  = nx*ny

    CALL GPartComm_Copy2Reg(this,vx,vxext)
    CALL GPartComm_Copy2Reg(this,vy,vyext)
    CALL GPartComm_Copy2Reg(this,vz,vzext)
    IF ( method .EQ. UNPACK_REP ) THEN
      DO k = 1, ngz  ! bottom extended zones
        km = k-1
        DO j=1,ny
          jm = j-1
          DO i=1,nx
            ! set top bcs:
            vx(i+(j-1)*nx+(k-1)*nxy) = vxext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy)
            vy(i+(j-1)*nx+(k-1)*nxy) = vyext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy)
            vz(i+(j-1)*nx+(k-1)*nxy) = vzext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy)

            ! set bottom bcs:
            vx(i+(j-1)*nx+(nz-ngz+k-1)*nxy) = vxext(i+ngp+(j+ngp-1)*nex+    (k-1)*nexy)
            vy(i+(j-1)*nx+(nz-ngz+k-1)*nxy) = vyext(i+ngp+(j+ngp-1)*nex+    (k-1)*nexy)
            vz(i+(j-1)*nx+(nz-ngz+k-1)*nxy) = vzext(i+ngp+(j+ngp-1)*nex+    (k-1)*nexy)

          ENDDO
        ENDDO
      ENDDO
    ELSE IF ( method .EQ. UNPACK_SUM ) THEN
      DO k = 1, ngz  ! bottom extended zones
        km = k-1
        DO j=1,ny
          jm = j-1
          DO i=1,nx
            ! set top bcs:
            vx(i+(j-1)*nx+(k-1)*nxy) = vx(i+(j-1)*nx+(k-1)*nxy) &
                              + vxext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy)
            vy(i+(j-1)*nx+(k-1)*nxy) = vy(i+(j-1)*nx+(k-1)*nxy) &
                              + vyext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy)
            vz(i+(j-1)*nx+(k-1)*nxy) = vz(i+(j-1)*nx+(k-1)*nxy) &
                              + vzext(i+ngp+(j+ngp-1)*nex+(nez+k-1)*nexy)

            ! set bottom bcs:
            vx(i+(j-1)*nx+(nz-ngz+k-1)*nxy) = vx(i+(j-1)*nx+(nz-ngz+k-1)*nxy) &
                              + vxext(i+ngp+(j+ngp-1)*nex+    (k-1)*nexy)
            vy(i+(j-1)*nx+(nz-ngz+k-1)*nxy) = vy(i+(j-1)*nx+(nz-ngz+k-1)*nxy) &
                              + vyext(i+ngp+(j+ngp-1)*nex+    (k-1)*nexy)
            vz(i+(j-1)*nx+(nz-ngz+k-1)*nxy) = vz(i+(j-1)*nx+(nz-ngz+k-1)*nxy) &
                              + vzext(i+ngp+(j+ngp-1)*nex+    (k-1)*nexy)

          ENDDO
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE GPartComm_LocalDataRetMF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_Copy2Reg(this,v,vext)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Copy2Reg
!  DESCRIPTION: Copy field from extended to regular grid,
!               without ghost zones.
!
!  ARGUMENTS  :
!    this     : 'this' class instance (IN)
!    v        : regular-grid field
!    vext     : extended-grid field
!    ldims    : local dims of v
!-----------------------------------------------------------------
!$  USE threads
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)                      :: this
    INTEGER                                             :: i,j,jm,k,km,ngp,ngz,nex,nexy
    INTEGER                                             :: nx,nxy,ny
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: v
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: vext

    ngz  = this%nzghost_
    ngp  = ngz * this%iextperp_
    nexy = (this%nd_(1)+2*ngp) * (this%nd_(2)+2*ngp)
    nex  = this%nd_(1)+2*ngp
    nx   = this%nd_(1)
    ny   = this%nd_(2)
    nxy  = nx*ny
!$omp parallel do if(ny*(this%kend_-this%ksta_+1).GE.nth) collapse(2) private(jm,i)
    DO km = 0,this%kend_-this%ksta_
!    DO k = 1,this%kend_-this%ksta_+1
!      km = k-1
      DO j=1,ny
        jm = j-1
        DO i=1,nx
          v(i+jm*nx+km*nxy) = vext(i+ngp+(jm+ngp)*nex+(km+ngz)*nexy)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE GPartComm_Copy2Reg
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!!!
END MODULE class_GPartComm
