!=================================================================
! GHOST GPartComm particles communication class. It handles 
!       two types of exchanges: the particles, and the 
!       velocity data used to update the particle positions,
!       (the 'ghost' zone data) and separate interfaces are 
!       provided  for each.
!       
!
! 2011 D. Rosenberg
!      ORNL: NCCS
!
! 15 Aug 2011: Initial version
!=================================================================
MODULE class_GPartComm
      USE mpivars
      USE fprecision
      IMPLICIT NONE
      
      TYPE, PUBLIC :: GPartComm
        PRIVATE
        ! Member data:
        INTEGER                                      :: maxparts_,nbuff_  ,nd_(3)   ,nzghost_
        INTEGER                                      :: nbsnd_   ,ntsnd_  ,nbrcv_   ,ntrcv_
        INTEGER                                      :: nprocs_  ,myrank_ ,comm_
        INTEGER                                      :: ntop_    ,nbot_   ,ierr_    ,istatus_
        INTEGER, ALLOCATABLE, DIMENSION(:,:)         :: ibsnd_   ,itsnd_  ,ibrcv_   ,itrcv_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: ibsh_    ,itsh_   ,ibrh_    ,itrh_ 
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: ibsndnz_ ,itsndnz_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: itop_    ,ibot_
        INTEGER, ALLOCATABLE, DIMENSION(:,:)         :: ibsnddst_,itsnddst_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: ibrcvnz_ ,itrcvnz_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION(:)     :: sbbuff_  ,stbuff_ ,rbbuff_  ,rtbuff_
      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: GPartComm_ctor
        FINAL            :: GPartComm_dtor
!!      PROCEDURE, PUBLIC:: GPartComm_dtor
        PROCEDURE,PUBLIC :: Init              => GPartComm_Init
        PROCEDURE,PUBLIC :: PartExchange      => GPartComm_PartExchange
        PROCEDURE,PUBLIC :: SlabDataExchange  => GPartComm_SlabDataExchangeMF
      END TYPE GPartComm

      PRIVATE :: GPartComm_Init              
      PRIVATE :: GPartComm_SlabDataExchangeMF, GPartComm_PartExchange    
      PRIVATE :: GPartComm_PackMF            , GPartComm_UnpackMF
      PRIVATE :: GPartComm_PackP             , GPartComm_UnpackP 
      PRIVATE :: GPartComm_LocalDataExch


! Methods:
  CONTAINS

  SUBROUTINE GPartComm_ctor(this, maxparts, nd, nzghost, comm)
!-----------------------------------------------------------------
!  Main explicit constructor
!  ARGUMENTS:
!    this    : 'this' class instance
!    nparts  : no. particles allowed on grid. 
!    nd(3)   : x- ,y- , and z- dimensions of data
!    nzghost : 'z' : no. slices of each slab required to 
!              build 'ghost' zones.  If there'sfewer slices on 
!              adjacent tasks, method will go next tasks to find 
!              the information to fill ghost zones.
!    comm    : MPI communicator
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT):: this
    INTEGER, INTENT(IN)           :: maxparts, nd(3), nzghost
    INTEGER, INTENT(IN)           :: comm

    this%maxparts_  = maxparts
    this%nd_        = nd
    this%nzghost_   = nzghost
    this%comm_      = comm;
    CALL MPI_COMM_SIZE(this%comm_,this%nprocs_,this%ierr_)
    CALL MPI_COMM_RANK(this%comm_,this%myrank_,this%ierr_)


    this%nbuff_     = MAX(maxparts*3*GP+1,3*nd(1)*nd(2)*nzghost*GP+nzghost+1) + 1
    ALLOCATE(this%sbbuff_ (this%nbuff_))
    ALLOCATE(this%stbuff_ (this%nbuff_))
    ALLOCATE(this%rbbuff_ (this%nbuff_))
    ALLOCATE(this%rtbuff_ (this%nbuff_))
    ALLOCATE(this%ibot_   (maxparts))
    ALLOCATE(this%itop_   (maxparts))

  END SUBROUTINE GPartComm_ctor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_dtor(this)
!-----------------------------------------------------------------
!  Main explicit destructor
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    TYPE(GPartComm),INTENT(INOUT)        :: this

    IF ( ALLOCATED(this%sbbuff_) ) DEALLOCATE(this%sbbuff_)
    IF ( ALLOCATED(this%stbuff_) ) DEALLOCATE(this%stbuff_)
    IF ( ALLOCATED(this%rbbuff_) ) DEALLOCATE(this%rbbuff_)
    IF ( ALLOCATED(this%rtbuff_) ) DEALLOCATE(this%rtbuff_)
    IF ( ALLOCATED  (this%itop_) ) DEALLOCATE(this%itop_)
    IF ( ALLOCATED  (this%ibot_) ) DEALLOCATE(this%ibot_)
  
    IF ( ALLOCATED(this%ibrcv_) ) DEALLOCATE(this%ibrcv_)
    IF ( ALLOCATED(this%itrcv_) ) DEALLOCATE(this%itrcv_)
    IF ( ALLOCATED(this%ibsnd_) ) DEALLOCATE(this%ibsnd_)
    IF ( ALLOCATED(this%itsnd_) ) DEALLOCATE(this%itsnd_)
    IF ( ALLOCATED (this%ibrh_) ) DEALLOCATE(this%ibrh_)
    IF ( ALLOCATED (this%itrh_) ) DEALLOCATE(this%itrh_)
    IF ( ALLOCATED (this%ibsh_) ) DEALLOCATE(this%ibsh_)
    IF ( ALLOCATED (this%itsh_) ) DEALLOCATE(this%itsh_)
    IF ( ALLOCATED (this%ibsnddst_) ) DEALLOCATE(this%ibsnddst_)
    IF ( ALLOCATED (this%itsnddst_) ) DEALLOCATE(this%itsnddst_)

  END SUBROUTINE GPartComm_dtor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_Init(this)
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
 

    IF ( ALLOCATED(this%ibrcv_) ) DEALLOCATE(this%ibrcv_)
    IF ( ALLOCATED(this%itrcv_) ) DEALLOCATE(this%itrcv_)
    IF ( ALLOCATED(this%ibsnd_) ) DEALLOCATE(this%ibsnd_)
    IF ( ALLOCATED(this%itsnd_) ) DEALLOCATE(this%itsnd_)
    IF ( ALLOCATED (this%ibrh_) ) DEALLOCATE(this%ibrh_)
    IF ( ALLOCATED (this%itrh_) ) DEALLOCATE(this%itrh_)
    IF ( ALLOCATED (this%ibsh_) ) DEALLOCATE(this%ibsh_)
    IF ( ALLOCATED (this%itsh_) ) DEALLOCATE(this%itsh_)
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
    ALLOCATE(this%ibsnddst_(nt,this%nzghost_+1))
    ALLOCATE(this%itsnddst_(nt,this%nzghost_+1))

    ! Initialize all task/neighbor  lists with -1:
    this%ibrcv_   =-1; this%itrcv_   =-1; this%ibsnd_  =-1; this%itsnd_ =-1;
    this%ibsnddst_=-1; this%itsnddst_=-1;

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

  SUBROUTINE GPartComm_LocalDataExch(this,vxext,vyext,vzext,vx,vy,vz)
!-----------------------------------------------------------------
!  METHOD     : LocalDataExch
!  DESCRIPTION: Does 'bdy exchange' of velocity component, when there's
!               only a single MPI task.
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
    INTEGER                                             :: i,j,k,ng,ney,nexy,nez,nxy
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)            :: vx,vy,vz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)            :: vxext,vyext,vzext

    ng   = this%nzghost_
    nxy  = this%nd_(1) * this%nd_(2)
    nexy = (this%nd_(1)+ng) * (this%nd_(2)+ng)
    ney  = this%nd_(2)+ng
    nez  = this%nd_(3)+ng

    DO k = 1, ng  ! bottom extended zones
      DO j=1,this%nd_(2)
        DO i=1,this%nd_(1)
          ! set bottom bcs:
          vxext(i+ng+(j+ng)*ney+k*nexy) = vx(i+j*this%nd_(1)+(k-ng-1)*nxy)
          vyext(i+ng+(j+ng)*ney+k*nexy) = vy(i+j*this%nd_(1)+(k-ng-1)*nxy)
          vzext(i+ng+(j+ng)*ney+k*nexy) = vz(i+j*this%nd_(1)+(k-ng-1)*nxy)

          ! set top bcs:
          vxext(i+ng+(j+ng)*ney+(nez+k)*nexy) = vx(i+j*this%nd_(1)+(k+1)*nxy)
          vyext(i+ng+(j+ng)*ney+(nez+k)*nexy) = vy(i+j*this%nd_(1)+(k+1)*nxy)
          vzext(i+ng+(j+ng)*ney+(nez+k)*nexy) = vz(i+j*this%nd_(1)+(k+1)*nxy)
        ENDDO
      ENDDO
    ENDDO
   
  END SUBROUTINE GPartComm_LocalDataExch
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPartComm_SlabDataExchangeMF(this,vxext,vyext,vzext,vx,vy,vz)
!-----------------------------------------------------------------
!  METHOD     : Exchange
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
    CALL GPartComm_UnpackMF(this,vxext,vyext,vzext,this%rbbuff_,'b')
    CALL GPartComm_UnpackMF(this,vxext,vyext,vzext,this%rtbuff_,'t')


  END SUBROUTINE GPartComm_SlabDataExchangeMF
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPartComm_PackMF(this,buff,vx,vy,vz,isnd,sdir)
!-----------------------------------------------------------------
!  METHOD     : PackPMP
!  DESCRIPTION: packs snd buffer with fields
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
    INTEGER                                 :: i,j,k,m,nt,nxy
    REAL(KIND=GP),INTENT  (OUT),DIMENSION(*):: buff
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*):: vx,vy,vz
    CHARACTER*(*),INTENT   (IN)             :: sdir


    IF ( sdir(1:1).NE.'b' .AND. sdir(1:1).NE.'B' &
    .AND.sdir(1:1).NE.'t' .AND. sdir(1:1).NE.'T' ) THEN
      WRITE(*,*) 'GPartComm_PackMF: Bad direction descriptor'
      STOP
    ENDIF

    nxy = this%nd_(2) * this%nd_(1);
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
        DO j = 1, this%nd_(2)
          DO i = 1, this%nd_(1)
            buff(nt) = vx(i+j*this%nd_(1)+k*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO
    
      DO m = 1, this%ibsndnz_(isnd)
        k = this%ibsnd_(isnd,m)
        DO j = 1, this%nd_(2)
          DO i = 1, this%nd_(1)
            buff(nt) = vy(i+j*this%nd_(1)+k*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

      DO m = 1,this%ibsndnz_(isnd)
        k = this%ibsnd_(isnd,m)
        DO j = 1, this%nd_(2)
          DO i = 1, this%nd_(1)
            buff(nt) = vz(i+j*this%nd_(1)+k*nxy)
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
        DO j = 1, this%nd_(2)
          DO i = 1, this%nd_(1)
            buff(nt) = vx(i+j*this%nd_(1)+k*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO
    
      DO m = 1,this%itsndnz_(isnd)
        k = this%itsnd_(isnd,m)
        DO j = 1, this%nd_(2)
          DO i = 1, this%nd_(1)
            buff(nt) = vy(i+j*this%nd_(1)+k*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO
    
      DO m = 1,this%itsndnz_(isnd)
        k = this%itsnd_(isnd,m)
        DO j = 1, this%nd_(2)
          DO i = 1, this%nd_(1)
            buff(nt) = vz(i+j*this%nd_(1)+k*nxy)
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO
    
    ENDIF
    

  END SUBROUTINE GPartComm_PackMF
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPartComm_UnpackMF(this,vxe,vye,vze,buff,sdir)
!-----------------------------------------------------------------
!  METHOD     : UnpackMF
!  DESCRIPTION: Unpacks recv buffer with into extended fields
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
    INTEGER                                 :: i,j,k,m,ng,nt,nxy
    REAL(KIND=GP),INTENT  (OUT),DIMENSION(*):: buff
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*):: vxe,vye,vze
    CHARACTER(len=1), INTENT  (IN)          :: sdir

    IF ( sdir(1:1).NE.'b' .AND. sdir(1:1).NE.'B' &
    .AND.sdir(1:1).NE.'t' .AND. sdir(1:1).NE.'T' ) THEN
      WRITE(*,*) 'GPartComm_UnpackMF: Bad direction descriptor'
      STOP
    ENDIF

    nxy = this%nd_(1)*this%nd_(2)
    ng  = this%nzghost_;
    IF      ( sdir(1:1) .EQ. 'b' .OR. sdir(1:1) .EQ. 'B' ) THEN
    ! Pack from bottom buffer:
      nt = 1
      DO m = 1,int(buff(1))
        k = int(buff(m+1))
        DO j = 1, this%nd_(2)
          DO i = 1, this%nd_(1)
            vxe(i+ng+(j+ng)*this%nd_(1)+k*nxy) = buff(nt) 
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO
    
      DO m = 1,int(buff(1))
        k = int(buff(m+1))
        DO j = 1, this%nd_(2)
          DO i = 1, this%nd_(1)
            vye(i+ng+(j+ng)*this%nd_(1)+k*nxy) = buff(nt) 
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

      DO m = 1,int(buff(1))
        k = int(buff(m+1))
        DO j = 1, this%nd_(2)
          DO i = 1, this%nd_(1)
            vze(i+ng+(j+ng)*this%nd_(1)+k*nxy) = buff(nt) 
            nt = nt + 1
          ENDDO
        ENDDO
      ENDDO

    ELSE ! Pack from top buffer:


    ENDIF

  END SUBROUTINE GPartComm_UnpackMF
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPartComm_PartExchange(this,px,py,pz,id,nparts,zmin,zmax)
!-----------------------------------------------------------------
!  METHOD     : PartExchange
!  DESCRIPTION: Carries out particle exchange. Particles will
!               be re-ordered after this call. 
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    px,py,pz: x,y,z positions of particles on global grid
!    id      : particle global ids
!    nparts  : number of particles in arrays id, px,py,pz
!    zmin/max: min/max z-dimensions of current MPI task
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER      ,INTENT(INOUT)                :: nparts
    INTEGER      ,INTENT(INOUT),DIMENSION(*)   :: id
    INTEGER                                    :: j,ibrank,itrank
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)   :: px,py,pz
    REAL(KIND=GP),INTENT   (IN)                :: zmin,zmax

    IF ( this%nprocs_ .EQ. 1 ) RETURN ! nothing to do

    itrank = mod(this%myrank_,this%nprocs_)
    ibrank = this%myrank_-1
    IF ( ibrank.LT.0 ) ibrank = nprocs-1

    ! Find pointers into particle lists for parts that must
    ! be sent to the top and bottom tasks:
    this%nbot_ = 0
    this%ntop_ = 0
    DO j = 0, nparts
      IF ( pz(j).LT.zmin .OR. pz(j).GT.zmax) THEN ! bottom
        this%ibot_(this%nbot_+1) = j
        this%nbot_ = this%nbot_ + 1
      ENDIF
      IF ( pz(j).GT.zmax .OR. pz(j).LT.zmin) THEN ! top
        this%itop_(this%ntop_+1) = j
        this%ntop_ = this%ntop_ + 1
      ENDIF
    ENDDO

    ! Post receives:
    CALL MPI_IRECV(this%rbbuff_,this%nbuff_,GC_REAL,ibrank, &
                   this%comm_,this%ibrh_(1),this%ierr_)
    CALL MPI_IRECV(this%rtbuff_,this%nbuff_,GC_REAL,itrank, &
                   this%comm_,this%itrh_(1),this%ierr_)

    !
    ! send data:
    CALL GPartComm_PackP(this,this%sbbuff_,this%nbuff_,px,py,pz,id,nparts,this%ibot_,this%nbot_)
    CALL MPI_ISEND(this%sbbuff_,this%nbuff_,GC_REAL,ibrank, &
                   this%comm_,this%itsh_(1),this%ierr_)
    CALL GPartComm_PackP(this,this%sbbuff_,this%nbuff_,px,py,pz,id,nparts,this%itop_,this%ntop_)
    CALL MPI_ISEND(this%stbuff_,this%nbuff_,GC_REAL,itrank, &
                   this%comm_,this%itsh_(1),this%ierr_)

    ! Concatenate partcle list to remove particles sent away:
    CALL GPartComm_Concat(this,px,py,pz,id,nparts,this%ibot_,&
                          this%nbot_,this%itop_,this%ntop_)

    CALL MPI_WAIT(this%ibrh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%ibrh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%ibsh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%itsh_(1),this%istatus_,this%ierr_)

    ! Update particle list:
    CALL GPartComm_UnpackP(this,px,py,pz,id,nparts,this%rbbuff_,this%nbuff_)
    CALL GPartComm_UnpackP(this,px,py,pz,id,nparts,this%rtbuff_,this%nbuff_)

  END SUBROUTINE GPartComm_PartExchange
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPartComm_PackP(this,buff,nbuff,px,py,pz,id,nparts,iind,nind)
!-----------------------------------------------------------------
!  METHOD     : PackP
!  DESCRIPTION: Packs send buffer with particles
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    buff    : buffer into which to pack particles for sends
!    nbuff   : max buffer length
!    px,py,pz: x,y,z positions of particles on global grid
!    id      : particle global ids
!    nparts  : number of particles in arrays id, px,py,pz
!    iind    : pointers into px,py,px,id particle arrays for
!              particles to pack
!    nind    : no. particles to pack
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER      ,INTENT(INOUT)                :: nbuff,nparts,nind
    INTEGER      ,INTENT(INOUT),DIMENSION(*)   :: id,iind
    INTEGER                                    :: j,nb
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)   :: buff
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)   :: px,py,pz

    buff(1) = nind
    nb = 1
    DO j = 1, nind
      buff(nb+1) = id(iind(j))
      buff(nb+2) = px(iind(j))
      buff(nb+3) = py(iind(j))
      buff(nb+4) = pz(iind(j))
      nb = nb + 4
    ENDDO

  END SUBROUTINE GPartComm_PackP
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPartComm_Concat(this,px,py,pz,id,nparts,ibind,nbind,itind,ntind)
!-----------------------------------------------------------------
!  METHOD     : PackP
!  DESCRIPTION: Removes particles at indices itind,ibind,and
!               concatenates the partles lists
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    buff    : buffer into which to pack particles for sends
!    px,py,pz: x,y,z positions of particles on global grid,updated
!    id      : particle global ids, updated
!    nparts  : number of particles in arrays id, px,py,pz, 
!              updated
!    ibind   : list of indices of parts sent to bottom task
!    nbind   : no. indices in ibind
!    itind   : list of indices of parts sent to top
!    ntind   : no. indices in itind
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER      ,INTENT(INOUT)                :: nparts
    INTEGER      ,INTENT   (IN)                :: nbind,ntind
    INTEGER      ,INTENT(INOUT),DIMENSION(*)   :: id
    INTEGER      ,INTENT   (IN),DIMENSION(*)   :: ibind,itind
    INTEGER                                    :: i,j,ngood
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)   :: px,py,pz

    DO j = 1, nbind
      id(ibind(j)) = -1
    ENDDO
    DO j = 1, nbind
      id(itind(j)) = -1
    ENDDO

    ngood = nparts - (nbind+ntind)
    j     = 1
    DO i = 1, ngood
      DO WHILE ( j.LE.nparts .AND. id(j).LT.0 )
        j = j + 1
      ENDDO
      IF ( j.LE.nparts .AND. j.NE.i ) THEN         
        id(i) = id(j); id(j) = -1
        px(i) = px(j)
        py(i) = py(j)
        pz(i) = pz(j)
      ENDIF

    ENDDO
    nparts = ngood


  END SUBROUTINE GPartComm_Concat
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPartComm_UnpackP(this,px,py,pz,id,nparts,buff,nbuff)
!-----------------------------------------------------------------
!  METHOD     : UnpackP
!  DESCRIPTION: Unpacks recv buffer with particles. Partlcles
!               will be added directly to the existing particle list.
!               
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    px,py,pz: x,y,z positions of particles on global grid, updated
!              with new particles
!    id      : particle global ids, updated with new particles
!    nparts  : new number of particles in arrays id, px,py,pz, updated
!              with new particles
!    buff    : buffer from which particle data is read
!    nbuff   : buffer length
!    
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars
    IMPLICIT NONE

    CLASS(GPartComm),INTENT(INOUT)             :: this
    INTEGER      ,INTENT(INOUT)                :: nparts
    INTEGER      ,INTENT(INOUT),DIMENSION(*)   :: id
    INTEGER      ,INTENT   (IN)                :: nbuff
    INTEGER                                    :: j,nb
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*)   :: px,py,pz
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*)   :: buff

    nb = 1
    DO j = 1, buff(1)
      nparts = nparts + 1
      id(nparts) = int(buff(nb+1))
      px(nparts) =     buff(nb+2)
      py(nparts) =     buff(nb+3)
      pz(nparts) =     buff(nb+4)
      nb = nb+4 
    ENDDO

  END SUBROUTINE GPartComm_UnpackP
!-----------------------------------------------------------------
!-----------------------------------------------------------------


END MODULE class_GPartComm
