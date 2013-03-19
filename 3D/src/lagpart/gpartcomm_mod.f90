!=================================================================
! GHOST GPartComm particles communication class. It handles 
!       two types of exchanges: the particles, and the 
!       velocity data used to update the particle positions,
!       (the 'ghost' zone data) and separate interfaces are 
!       provided  for each. The velocity data can be exchanged 
!       for _any_ other field variable as well, although the
!       'multi-field' interfaces may not be appropriate for it.
!       
!
! 2013 D. Rosenberg
!      ORNL: NCCS
!
! 15 Aug 2011: Initial version
!=================================================================
MODULE class_GPartComm
      USE mpivars
      USE fprecision
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      
      INTEGER,PARAMETER                              :: GPNULL=-1 ! particle NULL value
      TYPE, PUBLIC :: GPartComm
        PRIVATE
        ! Member data:
        INTEGER                                      :: intrfc_   ! if >=1 use multi-field interface
        INTEGER                                      :: maxparts_  ,nbuff_     ,nd_(3)   ,nzghost_
        INTEGER                                      :: nbsnd_     ,ntsnd_     ,nbrcv_   ,ntrcv_
        INTEGER                                      :: nprocs_    ,myrank_    ,comm_    
        INTEGER                                      :: csize_     ,nstrip_
        INTEGER                                      :: ntop_      ,nbot_      ,ierr_    ,istatus_
        INTEGER                                      :: btransinit_,iextperp_  ,nth_
        INTEGER, ALLOCATABLE, DIMENSION(:,:)         :: ibsnd_     ,itsnd_     ,ibrcv_   ,itrcv_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: ibsh_      ,itsh_      ,ibrh_    ,itrh_ 
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: igsh_      ,igrh_  
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: itypes_    ,ityper_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: ibsndnz_   ,itsndnz_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: itop_      ,ibot_
        INTEGER, ALLOCATABLE, DIMENSION(:,:)         :: ibsnddst_  ,itsnddst_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: ibrcvnz_   ,itrcvnz_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION(:)     :: sbbuff_    ,stbuff_    ,rbbuff_  ,rtbuff_
      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: GPartComm_ctor
        FINAL            :: GPartComm_dtor
        PROCEDURE,PUBLIC :: Init              => GPartComm_Init
        PROCEDURE,PUBLIC :: GetNumGhost       => GPartComm_GetNumGhost
        PROCEDURE,PUBLIC :: GTranspose        => GPartComm_Transpose
        PROCEDURE,PUBLIC :: VDBSynch          => GPartComm_VDBSynch
        PROCEDURE,PUBLIC :: SetCacheParam     => GPartComm_SetCacheParam
        PROCEDURE,PUBLIC :: PartExchangePDB   => GPartComm_PartExchangePDB
        PROCEDURE,PUBLIC :: PartExchangeV     => GPartComm_PartExchangeV
        PROCEDURE,PUBLIC :: SlabDataExchangeMF=> GPartComm_SlabDataExchangeMF
        PROCEDURE,PUBLIC :: SlabDataExchangeSF=> GPartComm_SlabDataExchangeSF
        PROCEDURE,PUBLIC :: ConcatPDB         => GPartComm_ConcatPDB
        PROCEDURE,PUBLIC :: ConcatV           => GPartComm_ConcatV

        GENERIC  ,PUBLIC :: PartExchange      => PartExchangePDB,PartExchangeV
!       GENERIC  ,PUBLIC :: SlabDataExchange  => SlabDataExchangeMF,SlabDataExchangeSF
        GENERIC  ,PUBLIC :: Concat            => ConcatPDB,ConcatV
      END TYPE GPartComm

      PRIVATE :: GPartComm_Init              
      PRIVATE :: GPartComm_SlabDataExchangeMF, GPartComm_SlabDataExchangeSF
      PRIVATE :: GPartComm_LocalDataExchMF   , GPartComm_LocalDataExchSF
      PRIVATE :: GPartComm_PartExchangePDB   , GPartComm_PartExchangeV 
      PRIVATE :: GPartComm_Transpose         , GPartComm_GetNumGhost
      PRIVATE :: GPartComm_PackMF            , GPartComm_UnpackMF
      PRIVATE :: GPartComm_PackSF            , GPartComm_UnpackSF
      PRIVATE :: GPartComm_PPackPDB          , GPartComm_PUnpackPDB
      PRIVATE :: GPartComm_PPackV            , GPartComm_PUnpackV 
      PRIVATE :: GPartComm_SetCacheParam


! Methods:
  CONTAINS

  SUBROUTINE GPartComm_ctor(this,intrface,maxparts,nd,nzghost,comm)
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
!    csize   : cache-size for local transposes
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT):: this
    INTEGER, INTENT(IN)           :: intrface,maxparts,nd(3),nzghost
    INTEGER, INTENT(IN)           :: comm
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
!$    this%nth_ = omp_get_max_threads()

    CALL MPI_COMM_SIZE(this%comm_,this%nprocs_,this%ierr_)
    CALL MPI_COMM_RANK(this%comm_,this%myrank_,this%ierr_)
    this%btransinit_ = .FALSE.

    IF ( this%intrfc_ .GE. 1 ) THEN
      this%nbuff_  = MAX(maxparts*4*(GP),3*nd(1)*nd(2)*nzghost*GP+nzghost+1) + 1
    ELSE
      this%nbuff_  = MAX(maxparts*4*(GP),nd(1)*nd(2)*nzghost*GP+nzghost+1) + 1
    ENDIF
    ALLOCATE(this%sbbuff_ (this%nbuff_))
    ALLOCATE(this%stbuff_ (this%nbuff_))
    ALLOCATE(this%rbbuff_ (this%nbuff_))
    ALLOCATE(this%rtbuff_ (this%nbuff_))
    ALLOCATE   (this%ibot_   (maxparts))
    ALLOCATE   (this%itop_   (maxparts))

  END SUBROUTINE GPartComm_ctor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_dtor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Main explicit destructor
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    TYPE(GPartComm),INTENT(INOUT)        :: this

    IF ( ALLOCATED   (this%sbbuff_) ) DEALLOCATE  (this%sbbuff_)
    IF ( ALLOCATED   (this%stbuff_) ) DEALLOCATE  (this%stbuff_)
    IF ( ALLOCATED   (this%rbbuff_) ) DEALLOCATE  (this%rbbuff_)
    IF ( ALLOCATED   (this%rtbuff_) ) DEALLOCATE  (this%rtbuff_)
    IF ( ALLOCATED     (this%itop_) ) DEALLOCATE    (this%itop_)
    IF ( ALLOCATED     (this%ibot_) ) DEALLOCATE    (this%ibot_)
  
    IF ( ALLOCATED    (this%ibrcv_) ) DEALLOCATE   (this%ibrcv_)
    IF ( ALLOCATED    (this%itrcv_) ) DEALLOCATE   (this%itrcv_)
    IF ( ALLOCATED    (this%ibsnd_) ) DEALLOCATE   (this%ibsnd_)
    IF ( ALLOCATED    (this%itsnd_) ) DEALLOCATE   (this%itsnd_)
    IF ( ALLOCATED     (this%ibrh_) ) DEALLOCATE    (this%ibrh_)
    IF ( ALLOCATED     (this%itrh_) ) DEALLOCATE    (this%itrh_)
    IF ( ALLOCATED     (this%ibsh_) ) DEALLOCATE    (this%ibsh_)
    IF ( ALLOCATED     (this%itsh_) ) DEALLOCATE    (this%itsh_)
    IF ( ALLOCATED     (this%igrh_) ) DEALLOCATE    (this%igrh_)
    IF ( ALLOCATED     (this%igsh_) ) DEALLOCATE    (this%igsh_)
    IF ( ALLOCATED   (this%itypes_) ) DEALLOCATE  (this%itypes_)
    IF ( ALLOCATED   (this%ityper_) ) DEALLOCATE  (this%ityper_)
    IF ( ALLOCATED (this%ibsnddst_) ) DEALLOCATE(this%ibsnddst_)
    IF ( ALLOCATED (this%itsnddst_) ) DEALLOCATE(this%itsnddst_)

    CALL GPartComm_Init(this)

  END SUBROUTINE GPartComm_dtor
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
    INTEGER                               :: i,j,k,ibrank,itrank,ksta,kend
    INTEGER                               :: kbsta,kbend,ktsta,ktend,n2p,nt

    ! Compute the no. sendto and recv from tasks there are:
 

    IF ( ALLOCATED    (this%ibrcv_) ) DEALLOCATE   (this%ibrcv_)
    IF ( ALLOCATED    (this%itrcv_) ) DEALLOCATE   (this%itrcv_)
    IF ( ALLOCATED    (this%ibsnd_) ) DEALLOCATE   (this%ibsnd_)
    IF ( ALLOCATED    (this%itsnd_) ) DEALLOCATE   (this%itsnd_)
    IF ( ALLOCATED     (this%ibrh_) ) DEALLOCATE    (this%ibrh_)
    IF ( ALLOCATED     (this%itrh_) ) DEALLOCATE    (this%itrh_)
    IF ( ALLOCATED     (this%ibsh_) ) DEALLOCATE    (this%ibsh_)
    IF ( ALLOCATED     (this%itsh_) ) DEALLOCATE    (this%itsh_)
    IF ( ALLOCATED     (this%igrh_) ) DEALLOCATE    (this%igrh_)
    IF ( ALLOCATED     (this%igsh_) ) DEALLOCATE    (this%igsh_)
    IF ( ALLOCATED   (this%itypes_) ) DEALLOCATE  (this%itypes_)
    IF ( ALLOCATED   (this%ityper_) ) DEALLOCATE  (this%ityper_)
    IF ( ALLOCATED (this%ibsnddst_) ) DEALLOCATE(this%ibsnddst_)
    IF ( ALLOCATED (this%itsnddst_) ) DEALLOCATE(this%itsnddst_)

    ! If there aren't enough 'slices' with nearest neighbors to 
    ! fill ghost zones, go to next furthest tasks, etc., to fill:
    this%nbsnd_ = 0
    this%ntsnd_ = 0
    this%nbrcv_ = 0
    this%ntrcv_ = 0
    n2p = this%nd_(3)/this%nprocs_
    nt  = this%nzghost_/n2p + 1 ! max no. tasks needed for ghost zones

    ALLOCATE(this%ibrcv_(nt,this%nzghost_+1))
    ALLOCATE(this%itrcv_(nt,this%nzghost_+1))
    ALLOCATE(this%ibsnd_(nt,this%nzghost_+1))
    ALLOCATE(this%itsnd_(nt,this%nzghost_+1))
    ALLOCATE(this%ibrh_(nt))
    ALLOCATE(this%itrh_(nt))
    ALLOCATE(this%ibsh_(nt))
    ALLOCATE(this%itsh_(nt))
    ALLOCATE(this%igrh_(0:this%nprocs_-1))
    ALLOCATE(this%igsh_(0:this%nprocs_-1))
    ALLOCATE(this%itypes_(0:this%nprocs_-1))
    ALLOCATE(this%ityper_(0:this%nprocs_-1))
    ALLOCATE(this%ibsnddst_(nt,this%nzghost_+1))
    ALLOCATE(this%itsnddst_(nt,this%nzghost_+1))

    ! Initialize all task/neighbor  lists with GPNULL:
    this%ibrcv_   =GPNULL; this%itrcv_   =GPNULL; this%ibsnd_  =GPNULL; this%itsnd_ =GPNULL;
    this%ibsnddst_=GPNULL; this%itsnddst_=GPNULL;

    ! Get global z-bounds on this rank:
    CALL range(1,this%nd_(3),this%nprocs_,this%myrank_,ksta,kend)

    j = 1
    ! Loop over possible tasks to find neighbor list info:
    DO WHILE ( j.LE.nt .AND. this%nbsnd_.LT.this%nzghost_ ) 
      itrank = mod(this%myrank_+j,this%nprocs_)
      ibrank = this%myrank_-j
      IF ( ibrank.lt.0 ) ibrank = this%nprocs_-j+1
      CALL range(1,this%nd_(3),this%nprocs_,ibrank,kbsta,kbend)
      CALL range(1,this%nd_(3),this%nprocs_,itrank,ktsta,ktend)
      this%ibrcv_(j,1) = ibrank    !bottom
      this%itrcv_(j,1) = itrank    !top

      k = 1
      DO WHILE ( k.LE.kend-ksta+1 .AND. this%nbsnd_.LE.this%nzghost_ ) 
        ! local z-indices to send to bottom & top tasks:
        this%ibsnd_   (j,k+1) = k               ! local z-index to be sent
        this%itsnd_   (j,k+1) = kend-ksta+1-k+1 ! local z-index to be sent
        this%ibsnddst_(j,k  ) = kend-ksta+1 + this%nzghost_ + k
        this%itsnddst_(j,k  ) = k
        this%nbsnd_ = this%nbsnd_ + 1
        this%ntsnd_ = this%ntsnd_ + 1
        k = k + 1
      ENDDO

      k = 1
      DO WHILE ( k.LE.kbend-kbsta+1 .AND. this%nbrcv_.LE.this%nzghost_ ) 
        ! local z-indices received from bottom task;
        ! where to place locally in extended array:
        ! (dims of extended array are 1:kend-ksta+1+2*nzghost)
        this%ibrcv_  (j,k+1) = this%nzghost_-k+1
        this%nbrcv_      = this%nbrcv_ + 1
        k = k + 1
      ENDDO

      k = 1
      DO WHILE ( k.LE.ktend-ktsta+1 .AND. this%ntrcv_.LE.this%nzghost_ ) 
        ! local z-indices received from top task;
        ! where to place locally in extended array:
        ! (dims of extended array are 1:kend-ksta+1+2*nzghost)
        this%itrcv_(j,k+1) = (kend-ksta+1)+this%nzghost_+k
        this%ntrcv_ = this%ntrcv_ + 1
        k = k + 1
      ENDDO
      j = j + 1
    ENDDO

    IF ( this%nbrcv_ .NE. this%nzghost_ &
    .OR. this%nbsnd_ .NE. this%nzghost_ & 
    .OR. this%ntrcv_ .NE. this%nzghost_ & 
    .OR. this%ntsnd_ .NE. this%nzghost_  ) THEN
       WRITE(*,*) 'GPartComm_Init: data incompatible with interpolation order' 
       STOP
    ENDIF

    ! Collect some data for easier sends, recvs:
    DO j=1,nt
      this%ibrcvnz_ = 0
      this%itrcvnz_ = 0
      this%ibsndnz_ = 0
      this%itsndnz_ = 0
      DO k=2,this%nzghost_
        IF ( this%ibsnd_(j,k).GE.0 ) this%ibsndnz_(j) = this%ibsndnz_(j) + 1
        IF ( this%itsnd_(j,k).GE.0 ) this%itsndnz_(j) = this%itsndnz_(j) + 1
        IF ( this%ibrcv_(j,k).GE.0 ) this%ibrcvnz_(j) = this%ibrcvnz_(j) + 1
        IF ( this%itrcv_(j,k).GE.0 ) this%itrcvnz_(j) = this%itrcvnz_(j) + 1
      ENDDO
    ENDDO


  END SUBROUTINE GPartComm_Init
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
!                        this%iextperp_=1 in contructor, but, again, these are
!                        bot filled here..
!    vx,vy,vz          : Eulerian velocity components on regular grid. Must
!                        be of size nd_ set in constructor
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)                      :: this
    INTEGER                                             :: i,j,k,ngp,ngz,ney,nexy,nez
    INTEGER                                             :: nx,nxy,ny,nz
    INTEGER                                             :: jm,km
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: vx,vy,vz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: vxext,vyext,vzext

    ngz  = this%nzghost_
    ngp  = ngz * this%iextperp_
    nexy = (this%nd_(1)+2*ngp) * (this%nd_(2)+2*ngp)
    ney  = this%nd_(2)+2*ngp
    nx   = this%nd_(1)
    ny   = this%nd_(2)
    nz   = this%nd_(3)
    nxy  = nx*ny

    DO k = 1, ngz  ! bottom extended zones
      km = k-1
      DO j=1,ny
        jm = j-1
        DO i=1,nx
          ! set bottom bcs:
          vxext(i+ngp+(j+ngp-1)*ney+    (k-1)*nexy) = vx(i+(j-1)*nx+(nz-ngz+k-2)*nxy)
          vyext(i+ngp+(j+ngp-1)*ney+    (k-1)*nexy) = vy(i+(j-1)*nx+(nz-ngz+k-2)*nxy)
          vzext(i+ngp+(j+ngp-1)*ney+    (k-1)*nexy) = vz(i+(j-1)*nx+(nz-ngz+k-2)*nxy)

          ! set top bcs:
          vxext(i+ngp+(j+ngp-1)*ney+(nez+k-1)*nexy) = vx(i+(j-1)*nx+k*nxy)
          vyext(i+ngp+(j+ngp-1)*ney+(nez+k-1)*nexy) = vy(i+(j-1)*nx+k*nxy)
          vzext(i+ngp+(j+ngp-1)*ney+(nez+k-1)*nexy) = vz(i+(j-1)*nx+k*nxy)

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
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)                      :: this
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: vx,vy,vz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: vxext,vyext,vzext

    INTEGER                                             :: itask,j,m
    IF ( this%intrfc_ .LT. 1 ) THEN
      WRITE(*,*) 'GPartComm_SlabDataExchangeMF: SF interface expected'
      STOP
    ENDIF

    IF ( this%nprocs_ .EQ. 1 ) THEN
      CALL GPartComm_LocalDataExch(this,vxext,vyext,vzext,vx,vy,vz)
      RETURN
    ENDIF

    ! post receives:
    DO m=1,this%nbrcv_  ! from bottom task:
      itask = this%ibrcv_(m,1)
      CALL MPI_IRECV(this%rbbuff_(this%ibrcv_(m,2)),this%nbuff_,GC_REAL,itask, &
                     this%comm_,this%ibrh_(m),this%ierr_) 
    ENDDO

    DO m=1,this%ntrcv_  ! from top task:
      itask = this%itrcv_(m,1)
      CALL MPI_IRECV(this%rtbuff_(this%itsnddst_(m,2)),this%nbuff_,GC_REAL,itask, &
                     this%comm_,this%itrh_(m),this%ierr_) 
    ENDDO


    !
    ! send data:
    DO m=1,this%nbsnd_  ! to bottom task:
      itask = this%ibsnd_(m,1)
      CALL GPartComm_PackMF(this,this%sbbuff_,vx,vy,vz,m,'b')
      CALL MPI_ISEND(this%sbbuff_,this%nbuff_,GC_REAL,this%itsnd_(m,1), &
                     this%comm_,this%itsh_(m),this%ierr_)  
    ENDDO
    DO m=1,this%ntsnd_  ! to top task:
      CALL GPartComm_PackMF(this,this%stbuff_,vx,vy,vz,m,'t')
      CALL MPI_ISEND(this%stbuff_,this%nbuff_,GC_REAL,this%itsnd_(m,1), &
                     this%comm_,this%itsh_(m),this%ierr_) 
    ENDDO

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

    ! Unpack received data:
    CALL GPartComm_UnpackMF(this,vxext,vyext,vzext,this%rbbuff_)
    CALL GPartComm_UnpackMF(this,vxext,vyext,vzext,this%rtbuff_)


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
    USE fprecision
    USE commtypes
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
        k = k-1
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
    USE fprecision
    USE commtypes
    USE mpivars
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)          :: this
    INTEGER                                 :: i,j,k,m,ngp,ngz,nt,nx,nxy,ny
    INTEGER                                 :: jm,km
    REAL(KIND=GP),INTENT  (OUT),DIMENSION(*):: buff
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*):: vxe,vye,vze

    nx  = this%nd_(1)
    ny  = this%nd_(2)
    nxy = nx*ny
    ngz = this%nzghost_;
    ngp = ngz*this%iextperp_

  ! Pack from either buffer:
    nt = 1
    DO m = 1,int(buff(1))
      k = int(buff(m+1))
      k = k-1
      DO j = 1, ny
        jm = j-1
        DO i = 1, nx
          vxe(i+ngp+(jm+ngp)*nx+km*nxy) = buff(nt)
          nt = nt + 1
        ENDDO
      ENDDO
    ENDDO

    DO m = 1,int(buff(1))
      k = int(buff(m+1))
      km = k-1
      DO j = 1, ny
        jm = j-1
        DO i = 1, nx
          vye(i+ngp+(jm+ngp)*nx+km*nxy) = buff(nt)
          nt = nt + 1
        ENDDO
      ENDDO
    ENDDO

    DO m = 1,int(buff(1))
      k = int(buff(m+1))
      km = k-1
      DO j = 1, ny
        jm = j-1
        DO i = 1, nx
          vze(i+ngp+(jm+ngp)*nx+km*nxy) = buff(nt)
          nt = nt + 1
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
!               only a single MPI task.. This is a single-field interface.
!  ARGUMENTS  :
!    this              : 'this' class instance (IN)
!    vxext,vyext,vzext : Eulerian velocity components returned on extended
!                        grid (that used to hold ghost zones). Only z-conditions
!                        are imposed; lateral periodicity is not handled here.
!    vx,vy,vz          : Eulerian velocity components on regular grid. Must
!                        be of size nd_ set in constructor
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)                      :: this
    INTEGER                                             :: i,j,k,ngp,ngz,ney,nexy,nez
    INTEGER                                             :: nx,nxy,ny,nz
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: v
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: vext

    ngz  = this%nzghost_
    ngp  = ngz * this%iextperp_
    nexy = (this%nd_(1)+2*ngp) * (this%nd_(2)+2*ngp)
    ney  = this%nd_(2)+2*ngp
    nx   = this%nd_(1)
    ny   = this%nd_(2)
    nz   = this%nd_(3)
    nxy  = nx*ny

    DO k = 1, ngz  ! bottom extended zones
      DO j=1,ny
        DO i=1,nx
          ! set bottom bcs:
          vext(i+ngp+(j+ngp-1)*ney+    (k-1)*nexy) = v(i+(j-1)*nx+(nz-ngz+k-2)*nxy)

          ! set top bcs:
          vext(i+ngp+(j+ngp-1)*ney+(nez+k-1)*nexy) = v(i+(j-1)*nx+k*nxy)
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
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)                      :: this
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: v
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: vext

    INTEGER                                             :: itask,j,m

    IF ( this%intrfc_ .GE. 1 ) THEN
      WRITE(*,*) 'GPartComm_SlabDataExchangeSF: MF interface expected'
      STOP
    ENDIF

    IF ( this%nprocs_ .EQ. 1 ) THEN
      CALL GPartComm_LocalDataExchSF(this,vext,v)
      RETURN
    ENDIF

    ! post receives:
    DO m=1,this%nbrcv_  ! from bottom task:
      itask = this%ibrcv_(m,1)
      CALL MPI_IRECV(this%rbbuff_(this%ibrcv_(m,2)),this%nbuff_,GC_REAL,itask, &
                     this%comm_,this%ibrh_(m),this%ierr_) 
    ENDDO

    DO m=1,this%ntrcv_  ! from top task:
      itask = this%itrcv_(m,1)
      CALL MPI_IRECV(this%rtbuff_(this%itsnddst_(m,2)),this%nbuff_,GC_REAL,itask, &
                     this%comm_,this%itrh_(m),this%ierr_) 
    ENDDO


    !
    ! send data:
    DO m=1,this%nbsnd_  ! to bottom task:
      itask = this%ibsnd_(m,1)
      CALL GPartComm_PackSF(this,this%sbbuff_,v,m,'b')
      CALL MPI_ISEND(this%sbbuff_,this%nbuff_,GC_REAL,this%itsnd_(m,1), &
                     this%comm_,this%itsh_(m),this%ierr_)  
    ENDDO
    DO m=1,this%ntsnd_  ! to top task:
      CALL GPartComm_PackSF(this,this%stbuff_,v,m,'t')
      CALL MPI_ISEND(this%stbuff_,this%nbuff_,GC_REAL,this%itsnd_(m,1), &
                     this%comm_,this%itsh_(m),this%ierr_) 
    ENDDO

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

    ! Unpack received data:
    CALL GPartComm_UnpackSF(this,vext,this%rbbuff_)
    CALL GPartComm_UnpackSF(this,vext,this%rtbuff_)


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
    USE fprecision
    USE commtypes
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
            buff(nt) = v(i+jm*nx+km*nxy)
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
            buff(nt) = v(i+jm*nx+km*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO
    
    ENDIF
    

  END SUBROUTINE GPartComm_PackSF
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPartComm_UnpackSF(this,vext,buff)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : UnpackSF
!  DESCRIPTION: Unpacks recv buffer with into extended (single) field
!               Messages received are self-referential, so contain info
!               on where to 'send' recvd data. So, there is no 't' or
!               'b' designation required for unpacking.
!  ARGUMENTS  :
!    this        : 'this' class instance (IN)
!    buff        : packed buffer (input) from which to store into 
!                  extended grid quantities.
!    vext        : Eulerian velocity component on extended grid
!                  in phys. space (IN)
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)          :: this
    INTEGER                                 :: i,j,k,m,ngp,ngz,nt,nx,nxy,ny
    INTEGER                                 :: jm,km
    REAL(KIND=GP),INTENT  (OUT),DIMENSION(*):: buff
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*):: vext

    nx   = this%nd_(1)
    ny   = this%nd_(2)
    nxy  = nx*ny
    ngz  = this%nzghost_
    ngp  = ngz*this%iextperp_

    ! Pack from either buffer:
    nt = 1
    DO m = 1,int(buff(1))
      k = int(buff(m+1))
      km = k-1
      DO j = 1, ny
        jm = j-1
        DO i = 1, nx
          vext(i+ngp+(jm+ngp)*nx+km*nxy) = buff(nt) 
          nt = nt + 1
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE GPartComm_UnpackSF
!-----------------------------------------------------------------
!-----------------------------------------------------------------


!!  SUBROUTINE GPartComm_PartExchangePDB(this,pdb,nparts,zmin,zmax)
!!!-----------------------------------------------------------------
!!!-----------------------------------------------------------------
!!!  METHOD     : PartExchangePDB
!!!  DESCRIPTION: Carries out particle exchange. Particles will
!!!               be re-ordered after this call. Uses PDB interface.
!!!               Note: For this call to work, the particle positions
!!!               must not be periodized on entry. In the same way,
!!!               zmin/zmax must also not be periodized.
!!!
!!!               Note that here, a particle on either zmin or zmax
!!!               is considered to be outside the interval defined
!!!               by zmin/zmax.
!!!  ARGUMENTS  :
!!!    this    : 'this' class instance (IN)
!!!    pdb     : part. d.b.
!!!    nparts  : number of particles in pdb
!!!    zmin/max: min/max z-dimensions of current MPI task
!!!-----------------------------------------------------------------
!!    USE fprecision
!!    USE commtypes
!!    USE mpivars
!!    USE pdbtypes
!!    IMPLICIT NONE
!!
!!    CLASS(GPartComm),INTENT(INOUT)             :: this
!!    INTEGER      ,INTENT(INOUT)                :: nparts
!!    INTEGER                                    :: j,ibrank,itrank
!!    TYPE(GPDBrec),INTENT(INOUT),DIMENSION(*)   :: pdb
!!    REAL(KIND=GP),INTENT   (IN)                :: zmin,zmax
!!
!!    IF ( this%nprocs_ .EQ. 1 ) RETURN ! nothing to do
!!
!!    itrank = mod(this%myrank_,this%nprocs_)
!!    ibrank = this%myrank_-1
!!    IF ( ibrank.LT.0 ) ibrank = nprocs-1
!!
!!    ! Find pointers into particle lists for parts that must
!!    ! be sent to the top and bottom tasks:
!!    this%nbot_ = 0
!!    this%ntop_ = 0
!!    DO j = 0, nparts
!!      IF ( pdb(j)%z.LE.zmin ) THEN ! bottom
!!        this%nbot_ = this%nbot_ + 1
!!        this%ibot_(this%nbot_) = j
!!      ENDIF
!!      IF ( pdb(j)%z.GE.zmax ) THEN ! top
!!        this%ntop_ = this%ntop_ + 1
!!        this%itop_(this%ntop_) = j
!!      ENDIF
!!    ENDDO
!!
!!    ! Post receives:
!!    CALL MPI_IRECV(this%rbbuff_,this%nbuff_,GC_REAL,ibrank, &
!!                   this%comm_,this%ibrh_(1),this%ierr_)
!!    CALL MPI_IRECV(this%rtbuff_,this%nbuff_,GC_REAL,itrank, &
!!                   this%comm_,this%itrh_(1),this%ierr_)
!!
!!    !
!!    ! send data:
!!    CALL GPartComm_PPack(this,this%sbbuff_,this%nbuff_,pdb,nparts,this%ibot_,this%nbot_)
!!    CALL MPI_ISEND(this%sbbuff_,this%nbuff_,GC_REAL,ibrank, &
!!                   this%comm_,this%itsh_(1),this%ierr_)
!!    CALL GPartComm_PPack(this,this%sbbuff_,this%nbuff_,pdb,nparts,this%itop_,this%ntop_)
!!    CALL MPI_ISEND(this%stbuff_,this%nbuff_,GC_REAL,itrank, &
!!                   this%comm_,this%itsh_(1),this%ierr_)
!!
!!    ! Concatenate partcle list to remove particles sent away:
!!    CALL GPartComm_ConcatPDB(this,pdb,nparts,this%ibot_,&
!!                          this%nbot_,this%itop_,this%ntop_)
!!
!!    CALL MPI_WAIT(this%ibrh_(1),this%istatus_,this%ierr_)
!!    CALL MPI_WAIT(this%ibrh_(1),this%istatus_,this%ierr_)
!!    CALL MPI_WAIT(this%ibsh_(1),this%istatus_,this%ierr_)
!!    CALL MPI_WAIT(this%itsh_(1),this%istatus_,this%ierr_)
!!
!!    ! Update particle list:
!!    CALL GPartComm_PUnpack(this,pdb,nparts,this%rbbuff_,this%nbuff_)
!!    CALL GPartComm_PUnpack(this,pdb,nparts,this%rtbuff_,this%nbuff_)
!!
!!  END SUBROUTINE GPartComm_PartExchangePDB
!!!-----------------------------------------------------------------
!!!-----------------------------------------------------------------


  SUBROUTINE GPartComm_PartExchangeV(this,id,px,py,pz,nparts,zmin,zmax)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PartExchangeV
!  DESCRIPTION: Carries out particle exchange. Particles will
!               be re-ordered after this call. Uses V interface.
!               Note: For this call to work, the particle positions
!               must _not_ be periodized on entry. In the same way,
!               zmin/zmax must also _not_ be periodized.
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
!    id      : array of particle ids
!    px,py,px: arrays containing x,y,z positions of particles
!    nparts  : number of particles in pdb
!    zmin/max: min/max z-dimensions of current MPI task
!    gext    : (3,2) real array containing global grid extents (start and
!              stop boundaries in each direction).
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars
    USE pdbtypes
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)               :: this
    INTEGER      ,INTENT(INOUT)                  :: nparts
    INTEGER                                      :: j,ibrank,itrank
    INTEGER      ,INTENT(INOUT),DIMENSION(nparts):: id
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nparts):: px,py,pz
    REAL(KIND=GP),INTENT   (IN)                  :: zmin,zmax

    IF ( this%nprocs_ .EQ. 1 ) RETURN ! nothing to do

    itrank = mod(this%myrank_+1,this%nprocs_)
    ibrank = this%myrank_-1
    IF ( ibrank.LT.0 ) ibrank = nprocs-1

    ! Find pointers into particle lists for parts that must
    ! be sent to the top and bottom tasks:
    this%nbot_ = 0
    this%ntop_ = 0
    DO j = 0, nparts
      IF ( pz(j).LE.zmin ) THEN ! bottom
        this%nbot_ = this%nbot_ + 1
        this%ibot_(this%nbot_) = j
      ENDIF
      IF ( pz(j).GE.zmax ) THEN ! top
        this%ntop_ = this%ntop_ + 1
        this%itop_(this%ntop_) = j
      ENDIF
    ENDDO

    ! Post receives:
    CALL MPI_IRECV(this%rbbuff_,this%nbuff_,GC_REAL,ibrank, &
                   this%comm_,this%ibrh_(1),this%ierr_)
    CALL MPI_IRECV(this%rtbuff_,this%nbuff_,GC_REAL,itrank, &
                   this%comm_,this%itrh_(1),this%ierr_)

    !
    ! send data:
    CALL GPartComm_PPackV(this,this%sbbuff_,this%nbuff_,id,px,py,pz,nparts,this%ibot_,this%nbot_)
    CALL MPI_ISEND(this%sbbuff_,this%nbuff_,GC_REAL,ibrank, &
                   this%comm_,this%itsh_(1),this%ierr_)
    CALL GPartComm_PPackV(this,this%sbbuff_,this%nbuff_,id,px,py,pz,nparts,this%itop_,this%ntop_)
    CALL MPI_ISEND(this%stbuff_,this%nbuff_,GC_REAL,itrank, &
                   this%comm_,this%itsh_(1),this%ierr_)

    ! Concatenate partcle list to remove particles sent away:
    CALL GPartComm_ConcatV(this,id,px,py,pz,nparts,this%ibot_,&
                          this%nbot_,this%itop_,this%ntop_)

    CALL MPI_WAIT(this%ibrh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%ibrh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%ibsh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%itsh_(1),this%istatus_,this%ierr_)

    ! Update particle list:
    CALL GPartComm_PUnpackV(this,id,px,py,pz,nparts,this%rbbuff_,this%nbuff_)
    CALL GPartComm_PUnpackV(this,id,px,py,pz,nparts,this%rtbuff_,this%nbuff_)

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
    USE fprecision
    USE commtypes
    USE mpivars
    USE pdbtypes
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER      ,INTENT(INOUT)                :: nbuff,nparts,nind
    INTEGER      ,INTENT(INOUT),DIMENSION(*)   :: iind
    INTEGER      ,INTENT   (IN),DIMENSION(*)   :: id
    INTEGER                                    :: j,nb
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)   :: buff
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)   :: px,py,pz

    buff(1) = nind
    nb = 1
    DO j = 1, nind
      buff(nb+1) = id(iind(j))
      buff(nb+2) = px(iind(j))
      buff(nb+3) = py(iind(j))
      buff(nb+4) = pz(iind(j))
      nb = nb + 4
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
    USE fprecision
    USE commtypes
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)                :: this
    INTEGER      ,INTENT(INOUT)                   :: nparts
    INTEGER      ,INTENT   (IN)                   :: nbind,ntind
    INTEGER      ,INTENT   (IN),DIMENSION(nparts) :: ibind,itind
    INTEGER      ,INTENT(INOUT),DIMENSION(nparts) :: id
    INTEGER                                       :: i,j,ngood
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nparts) :: px,py,pz

    DO j = 1, nbind
      id(ibind(j)) = GPNULL
    ENDDO
    DO j = 1, nbind
      id(itind(j)) = GPNULL
    ENDDO

    ngood = nparts - (nbind+ntind)
    j     = 1
    DO i = 1, ngood
      DO WHILE ( j.LE.nparts .AND. id(j).EQ.GPNULL )
        j = j + 1
      ENDDO
      IF ( j.LE.nparts .AND. j.NE.i ) THEN         
        id(i) = id(j); id(j) = GPNULL
        px(i) = px(j)
        py(i) = py(j)
        pz(i) = pz(j)
      ENDIF

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
    USE fprecision
    USE commtypes
    USE mpivars
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)                :: this
    INTEGER      ,INTENT(INOUT)                   :: nparts
    INTEGER      ,INTENT   (IN)                   :: nbuff
    INTEGER      ,INTENT(INOUT),DIMENSION(nparts) :: id
    INTEGER                                       :: j,nb
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nparts) :: px,py,pz
    REAL(KIND=GP),INTENT   (IN),DIMENSION(nparts) :: buff

    nb = 1
    DO j = 1, buff(1)
      nparts = nparts + 1
      id(nparts) = int(buff(nb+1))
      px(nparts) =      buff(nb+2)
      py(nparts) =      buff(nb+3)
      pz(nparts) =      buff(nb+4)
      nb = nb+4 
    ENDDO

  END SUBROUTINE GPartComm_PUnpackV
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPartComm_PartExchangePDB(this,pdb,nparts,zmin,zmax)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : PartExchangePDB
!  DESCRIPTION: Carries out particle exchange. Particles will
!               be re-ordered after this call. Uses PDB  interface.
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
    USE fprecision
    USE commtypes
    USE mpivars
    USE pdbtypes
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER      ,INTENT(INOUT)                :: nparts
    INTEGER                                    :: j,ibrank,itrank
    TYPE(GPDBrec),INTENT(INOUT),DIMENSION(*)   :: pdb
    REAL(KIND=GP),INTENT   (IN)                :: zmin,zmax

    IF ( this%nprocs_ .EQ. 1 ) RETURN ! nothing to do

    itrank = mod(this%myrank_,this%nprocs_)
    ibrank = this%myrank_-1
    IF ( ibrank.LT.0 ) ibrank = nprocs-1

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
    CALL MPI_IRECV(this%rbbuff_,this%nbuff_,GC_REAL,ibrank, &
                   this%comm_,this%ibrh_(1),this%ierr_)
    CALL MPI_IRECV(this%rtbuff_,this%nbuff_,GC_REAL,itrank, &
                   this%comm_,this%itrh_(1),this%ierr_)

    !
    ! send data:
    CALL GPartComm_PPack(this,this%sbbuff_,this%nbuff_,pdb,nparts,this%ibot_,this%nbot_)
    CALL MPI_ISEND(this%sbbuff_,this%nbuff_,GC_REAL,ibrank, &
                   this%comm_,this%itsh_(1),this%ierr_)
    CALL GPartComm_PPack(this,this%sbbuff_,this%nbuff_,pdb,nparts,this%itop_,this%ntop_)
    CALL MPI_ISEND(this%stbuff_,this%nbuff_,GC_REAL,itrank, &
                   this%comm_,this%itsh_(1),this%ierr_)

    ! Concatenate partcle list to remove particles sent away:
    CALL GPartComm_ConcatPDB(this,pdb,nparts,this%ibot_,&
                          this%nbot_,this%itop_,this%ntop_)

    CALL MPI_WAIT(this%ibrh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%ibrh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%ibsh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%itsh_(1),this%istatus_,this%ierr_)

    ! Update particle list:
    CALL GPartComm_PUnpack(this,pdb,nparts,this%rbbuff_,this%nbuff_)
    CALL GPartComm_PUnpack(this,pdb,nparts,this%rtbuff_,this%nbuff_)

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
    USE fprecision
    USE commtypes
    USE mpivars
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
    USE fprecision
    USE commtypes
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
    USE fprecision
    USE commtypes
    USE mpivars
    USE pdbtypes
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER      ,INTENT(INOUT)                :: nparts
    INTEGER      ,INTENT   (IN)                :: nbuff
    INTEGER                                    :: j,nb
    TYPE(GPDBrec),INTENT(INOUT),DIMENSION(*)   :: pdb
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)   :: buff

    nb = 1
    DO j = 1, buff(1)
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


  SUBROUTINE GPartComm_Transpose(this,ofield,ifield,idims,rank,tmp)
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
!    ifield  : input field that is xy complete
!    idims   : local dimensions of ifield. 
!    rank    : rank of field (how many 'idims' array elements)
!    tmp     : real field of size required to hold field and its
!              transpose locally. Must be the same size on all 
!              MPI tasks
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars
    USE gtimer
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    REAL(KIND=GP),INTENT  (OUT),DIMENSION(*)   :: ofield
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)   :: ifield,tmp
    INTEGER      ,INTENT   (IN),DIMENSION(*)   :: idims
    INTEGER      ,INTENT   (IN)                :: rank
    INTEGER                                    :: i,ii,j,jj,k,kk
    INTEGER                                    :: igetfrom,iproc,irank,isendto,istrip
    INTEGER                                    :: nx,ny,nz,nxy,nzy


    IF ( .NOT.this%btransinit_ ) THEN
      IF ( rank.EQ.2 ) THEN
        CALL GPartComm_InitTrans2D(this)
      ENDIF
      IF ( rank.EQ.3 ) THEN
        CALL GPartComm_InitTrans3D(this)
      ENDIF
    ENDIF 

    ! NOTE: rank is transpose problem rank; irank is MPI rank...
    nx = idims(1)
    ny = idims(2)
    IF ( rank .GT. 2 ) nz = idims(3)
    nxy = nx*ny
    nzy = nz*ny

    DO iproc = 0, nprocs-1, this%nstrip_
       DO istrip=0, this%nstrip_-1
          irank = iproc + istrip

          isendto = this%myrank_ + irank
          if ( isendto .ge. nprocs ) isendto = isendto - nprocs

          igetfrom = this%myrank_- irank
          if ( igetfrom .lt. 0 ) igetfrom = igetfrom + nprocs
          CALL MPI_IRECV(tmp,1,this%ityper_(igetfrom),igetfrom,      &
                        1,this%comm_,this%igrh_(irank),this%ierr_)

          CALL MPI_ISEND(ifield,1,this%itypes_(isendto),isendto, &
                        1,this%comm_,this%igsh_(irank),this%ierr_)
       ENDDO

       DO istrip=0, this%nstrip_-1
          irank = iproc + istrip
          CALL MPI_WAIT(this%igsh_(irank),this%istatus_,this%ierr_)
          CALL MPI_WAIT(this%igrh_(irank),this%istatus_,this%ierr_)
       ENDDO
    ENDDO

    IF ( rank .EQ. 3 ) THEN

!$omp parallel do if ((idims(3)-1)/this%csize_.ge.this%nth_) private (jj,kk,i,j,k)
    DO ii = 1,nz,this%csize_
!$omp parallel do if ((idims(3)-1)/this%csize_.lt.this%nth_) private (kk,i,j,k)
       DO jj = 1,ny,this%csize_
          DO kk = 1,nx,this%csize_

             DO i = ii,min(nz,ii+this%csize_-1)
               DO j = jj,min(ny,jj+this%csize_-1)
                 DO k = kk,min(nx,kk+this%csize_-1)
                    ofield(k+(j-1)*nz+(i-1)*nzy) = tmp(i+(j-1)*nx+(k-1)*nxy)
                 END DO
               END DO
             END DO

          END DO
       END DO
    END DO

    ELSE


       DO ii = 1,ny,this%csize_
         DO jj = 1,nx,this%csize_
            DO i = ii,min(ny,ii+this%csize_-1)
              DO j = jj,min(nx,jj+this%csize_-1)
                 ofield(j+(i-1)*ny) = tmp(i+(j-1)*nx)
              END DO
            END DO
         END DO
      END DO

    ENDIF


  END SUBROUTINE GPartComm_Transpose
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
    USE commtypes
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER                                    :: ista,iend
    INTEGER                                    :: jsta,jend
    INTEGER                                    :: irank,jrank
    INTEGER                                    :: itemp1,itemp2
    
    CALL range(1,this%nd_(2),this%nprocs_,this%myrank_,jsta,jend)
    DO irank = 0,this%nprocs_-1
       CALL range(1,this%nd_(1),this%nprocs_,irank,ista,iend)
       CALL block2d(1,this%nd_(1),jsta,ista,iend,jsta,jend, &
                    GC_REAL,itemp1)
       this%itypes_(irank) = itemp1
    END DO
    CALL range(1,this%nd_(1),this%nprocs_,this%myrank_,ista,iend)
    DO jrank = 0,this%nprocs_-1
       CALL range(1,this%nd_(2),this%nprocs_,jrank,jsta,jend)
       CALL block2d(ista,iend,1,ista,iend,jsta,jend,  &
                   GC_REAL,itemp2)
       this%ityper_(jrank) = itemp2
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
    USE commtypes
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
       this%itypes_(irank) = itemp1
    END DO
    CALL range(1,this%nd_(1),this%nprocs_,this%myrank_,ista,iend)
    DO krank = 0,this%nprocs_-1
       CALL range(1,this%nd_(3),this%nprocs_,krank,ksta,kend)
       CALL block3d(ista,iend,1,this%nd_(2),1,ista,iend,1,this%nd_(2), &
                   ksta,kend,GC_REAL,itemp2)
       this%ityper_(krank) = itemp2
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
    USE commtypes
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)                  :: this
    INTEGER      ,INTENT   (IN),DIMENSION(*)        :: id
    INTEGER      ,INTENT   (IN)                     :: nl
    INTEGER      ,INTENT   (IN)                     :: ngvdb
    INTEGER                                         :: i,j
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)        :: lx,ly,lz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(3,ngvdb)  :: gvdb,ptmp


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
    CALL MPI_ALLREDUCE(gvdb,ptmp,3*ngvdb,GC_REAL,   &
                       MPI_SUM,this%comm_,this%ierr_)


  END SUBROUTINE GPartComm_VDBSynch
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

END MODULE class_GPartComm
