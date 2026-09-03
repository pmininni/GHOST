!=================================================================
! GPartComm: communication for the particle classes: exchange of
! the ghost z-planes of the field the particles interpolate from,
! exchange of the particles that leave the slab of a task,
! synchronization of the global particle database (VDB exchange),
! and the gather of the particle records on task 0 used by the
! non-collective I/O.
!
! In offload builds the field buffers, the particle send/receive
! buffers and the index lists have device copies (galloc); the
! kernels that fill or read them run on the device while
! gdev_active is set, and MPI receives the device addresses of the
! buffers (use_device_addr). Particle records travel as two
! contiguous messages, the ids and the three coordinates, which
! GPU-aware MPI handles without derived datatypes. The selection
! of particles (leaving through the bottom or the top of the slab,
! holes left by the departed particles) uses gpselect, so the
! local order of the particles is the same for any number of
! threads or teams.
!
! The routines used only at initialization or for the I/O
! (IdentifyTaskV, PartScatterV, VDBSynch_t0, LagSynch_t0) run on
! the host copies and must be called with gdev_active unset.
!=================================================================
MODULE class_GPartComm
      USE fprecision
      USE commtypes
      USE gtimer
      USE gmem
      USE gdevice, ONLY: gdev_active
      USE gpselect
      IMPLICIT NONE
      INTEGER,PARAMETER,PUBLIC :: GPNULL=-1          ! particle NULL value
      INTEGER,PARAMETER,PUBLIC :: GPCOMM_INTRFC_SF=0 ! single-field interface
      INTEGER,PARAMETER,PUBLIC :: GPEXCH_INIT = 0    ! starts particle exchange
      INTEGER,PARAMETER,PUBLIC :: GPEXCH_UPDT = 1    ! continues part. exchange
      INTEGER,PARAMETER,PUBLIC :: GPEXCH_END  = 2    ! finishes particle exchange
      INTEGER,PARAMETER,PUBLIC :: GPEXCH_UNIQ = 3    ! exchange only positions
      PRIVATE

      TYPE, PUBLIC :: GPartComm
        PRIVATE
        INTEGER                                      :: maxparts_  ,nbuff_     ,nd_(3)   ,nzghost_
        INTEGER                                      :: nbsnd_     ,ntsnd_     ,nbrcv_   ,ntrcv_
        INTEGER                                      :: nprocs_    ,myrank_
        INTEGER                                      :: ntop_      ,nbot_      ,ierr_
        INTEGER                                      :: iextperp_  ,ksta_      ,kend_    ,nth_
        INTEGER                                      :: hcomm_
        ! Transposes of a real field between the slab (nx,ny,kl) and the
        ! z-complete (nz,ny,il) layouts: x and z ranges of every task,
        ! offsets and counts of the contiguous messages, device buffers
        LOGICAL                                      :: btransinit_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: txsta_     ,txend_     ,tzsta_   ,tzend_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: tso_       ,tsc_       ,tro_     ,trc_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION(:)     :: gtsbuf_    ,gtrbuf_
        ! Tables of the ghost-plane exchange (host)
        INTEGER, ALLOCATABLE, DIMENSION(:,:)         :: ibsnd_     ,itsnd_
        INTEGER, ALLOCATABLE, DIMENSION(:,:)         :: ibsnddst_  ,itsnddst_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: ibrcv_     ,itrcv_     ,nbbrcv_  ,ntbrcv_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: ibsndnz_   ,itsndnz_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: ibrcvnz_   ,itrcvnz_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: ibsndp_    ,ibrcvp_    ,itsndp_  ,itrcvp_
        ! Index lists and buffers (device copies in offload builds)
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: itop_      ,ibot_      ,oldid_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: iflag_     ,ihole_     ,isurv_
        INTEGER, ALLOCATABLE, DIMENSION  (:)         :: sbid_      ,stid_      ,rbid_    ,rtid_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:)   :: sbpr_      ,stpr_      ,rbpr_    ,rtpr_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:)   :: sbbuff_    ,stbuff_
        REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:)   :: rbbuff_    ,rtbuff_
        TYPE(MPI_Comm)                               :: comm_
        TYPE(MPI_Status)                             :: istatus_
        TYPE(MPI_Request), DIMENSION(:), ALLOCATABLE :: ibsh_      ,ibrh_      ,itsh_    ,itrh_
      CONTAINS
        PROCEDURE,PUBLIC :: GPartComm_ctor
        FINAL            :: GPartComm_dtor
        PROCEDURE,PUBLIC :: Init              => GPartComm_Init
        PROCEDURE,PUBLIC :: GetNumGhost       => GPartComm_GetNumGhost
        PROCEDURE,PUBLIC :: GTranspose        => GPartComm_Transpose
        PROCEDURE,PUBLIC :: GITranspose       => GPartComm_ITranspose
        PROCEDURE,PUBLIC :: VDBSynch          => GPartComm_VDBSynch
        PROCEDURE,PUBLIC :: VDBSynch_t0       => GPartComm_VDBSynch_t0
        PROCEDURE,PUBLIC :: LagSynch_t0       => GPartComm_LagSynch_t0
        PROCEDURE,PUBLIC :: IdentifyTaskV     => GPartComm_IdentifyTaskV
        PROCEDURE,PUBLIC :: PartScatterV      => GPartComm_PartScatterV
        PROCEDURE,PUBLIC :: PartExchangeV     => GPartComm_PartExchangeV
        PROCEDURE,PUBLIC :: IdentifyExchV     => GPartComm_IdentifyExchV
        PROCEDURE,PUBLIC :: SlabDataExchangeSF=> GPartComm_SlabDataExchangeSF
        PROCEDURE,PUBLIC :: ResizeArrays      => GPartComm_ResizeArrays
      END TYPE GPartComm

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
    TYPE(mpi_comm), INTENT(IN)    :: comm
    INTEGER, INTENT(IN)           :: hcomm
!$  INTEGER, EXTERNAL             :: omp_get_max_threads

    IF ( intrface .NE. GPCOMM_INTRFC_SF ) THEN
      WRITE(*,*) 'GPartComm_ctor: only the single-field interface is supported'
      STOP
    ENDIF
    this%maxparts_  = maxparts
    this%nd_        = nd
    this%nzghost_   = nzghost
    this%comm_      = comm;
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

    this%nbuff_  = nd(1)*nd(2)*nzghost+nzghost+2 ! ghost planes + header

    CALL GPartComm_Init(this)

  END SUBROUTINE GPartComm_ctor


!-----------------------------------------------------------------
! Allocation of the particle-sized arrays (index lists, flags and
! send/receive buffers); called by Init
!-----------------------------------------------------------------
  SUBROUTINE GPartComm_AllocParts(this,np)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT) :: this
    INTEGER         ,INTENT(IN)    :: np
    CALL galloc(this%ibot_ ,np)
    CALL galloc(this%itop_ ,np)
    CALL galloc(this%oldid_,np)
    CALL galloc(this%iflag_,np)
    CALL galloc(this%ihole_,np)
    CALL galloc(this%isurv_,np)
    CALL galloc(this%sbid_ ,np)
    CALL galloc(this%stid_ ,np)
    CALL galloc(this%rbid_ ,np)
    CALL galloc(this%rtid_ ,np)
    CALL galloc(this%sbpr_ ,3,np)
    CALL galloc(this%stpr_ ,3,np)
    CALL galloc(this%rbpr_ ,3,np)
    CALL galloc(this%rtpr_ ,3,np)
  END SUBROUTINE GPartComm_AllocParts


  SUBROUTINE GPartComm_dtor(this)
    IMPLICIT NONE
    TYPE(GPartComm),INTENT(INOUT)        :: this
    CALL GPartComm_DoDealloc(this)
  END SUBROUTINE GPartComm_dtor


  SUBROUTINE GPartComm_DoDealloc(this)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)        :: this
    CALL gfree(this%sbbuff_); CALL gfree(this%stbuff_)
    CALL gfree(this%rbbuff_); CALL gfree(this%rtbuff_)
    CALL gfree(this%ibot_)  ; CALL gfree(this%itop_)  ; CALL gfree(this%oldid_)
    CALL gfree(this%iflag_) ; CALL gfree(this%ihole_) ; CALL gfree(this%isurv_)
    CALL gfree(this%sbid_)  ; CALL gfree(this%stid_)
    CALL gfree(this%rbid_)  ; CALL gfree(this%rtid_)
    CALL gfree(this%sbpr_)  ; CALL gfree(this%stpr_)
    CALL gfree(this%rbpr_)  ; CALL gfree(this%rtpr_)
    IF ( ALLOCATED(this%ibrcvp_)   ) DEALLOCATE(this%ibrcvp_)
    IF ( ALLOCATED(this%ibsndp_)   ) DEALLOCATE(this%ibsndp_)
    IF ( ALLOCATED(this%itrcvp_)   ) DEALLOCATE(this%itrcvp_)
    IF ( ALLOCATED(this%itsndp_)   ) DEALLOCATE(this%itsndp_)
    IF ( ALLOCATED(this%ibrcv_)    ) DEALLOCATE(this%ibrcv_)
    IF ( ALLOCATED(this%itrcv_)    ) DEALLOCATE(this%itrcv_)
    IF ( ALLOCATED(this%nbbrcv_)   ) DEALLOCATE(this%nbbrcv_)
    IF ( ALLOCATED(this%ntbrcv_)   ) DEALLOCATE(this%ntbrcv_)
    IF ( ALLOCATED(this%ibsnd_)    ) DEALLOCATE(this%ibsnd_)
    IF ( ALLOCATED(this%itsnd_)    ) DEALLOCATE(this%itsnd_)
    IF ( ALLOCATED(this%ibsnddst_) ) DEALLOCATE(this%ibsnddst_)
    IF ( ALLOCATED(this%itsnddst_) ) DEALLOCATE(this%itsnddst_)
    IF ( ALLOCATED(this%ibsndnz_)  ) DEALLOCATE(this%ibsndnz_)
    IF ( ALLOCATED(this%itsndnz_)  ) DEALLOCATE(this%itsndnz_)
    IF ( ALLOCATED(this%ibrcvnz_)  ) DEALLOCATE(this%ibrcvnz_)
    IF ( ALLOCATED(this%itrcvnz_)  ) DEALLOCATE(this%itrcvnz_)
    IF ( ALLOCATED(this%ibrh_)     ) DEALLOCATE(this%ibrh_)
    IF ( ALLOCATED(this%itrh_)     ) DEALLOCATE(this%itrh_)
    IF ( ALLOCATED(this%ibsh_)     ) DEALLOCATE(this%ibsh_)
    IF ( ALLOCATED(this%itsh_)     ) DEALLOCATE(this%itsh_)
    IF ( ALLOCATED(this%txsta_) ) DEALLOCATE(this%txsta_,this%txend_,this%tzsta_,this%tzend_)
    IF ( ALLOCATED(this%tso_)   ) DEALLOCATE(this%tso_,this%tsc_,this%tro_,this%trc_)
    CALL gfree(this%gtsbuf_); CALL gfree(this%gtrbuf_)
    this%btransinit_ = .FALSE.
  END SUBROUTINE GPartComm_DoDealloc


!=================================================================
! Ghost-plane exchange of a single field
!=================================================================

!-----------------------------------------------------------------
!  METHOD     : SlabDataExchangeSF
!  DESCRIPTION: Fills the extended field vext(nx,ny,nzl+2*nzghost)
!               with the local slab v(nx,ny,nzl) and with the
!               nzghost planes below and above it, received from
!               the neighbor tasks (periodic in z). Runs on the
!               device while gdev_active is set: the planes are
!               packed into device buffers, sent with GPU-aware
!               MPI and unpacked on the device. The messages are
!               self-describing: buff(1) holds the number of
!               planes, buff(2:npl+1) their destination indices in
!               the extended grid, and the planes follow.
!-----------------------------------------------------------------
  SUBROUTINE GPartComm_SlabDataExchangeSF(this,vext,v)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)           :: this
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*) :: v
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*) :: vext
    INTEGER                                  :: j,nx,ny,nzl,ngz,ngp,nex,ney,nez,nt

    nx  = this%nd_(1)
    ny  = this%nd_(2)
    nzl = this%kend_-this%ksta_+1
    ngz = this%nzghost_
    ngp = ngz*this%iextperp_
    nex = nx+2*ngp
    ney = ny+2*ngp
    nez = nzl+2*ngz
    nt  = SIZE(this%sbbuff_,2)

    CALL gpc_copy2ext(nx,ny,nzl,nex,ney,nez,ngp,ngz,v,vext)
    IF ( this%nprocs_ .EQ. 1 ) THEN
      CALL gpc_localexch(nx,ny,nzl,nex,ney,nez,ngp,ngz,v,vext)
      RETURN
    ENDIF

    DO j=1,this%nbsnd_  ! to bottom task:
      CALL gpc_packsf(nx,ny,nzl,this%ibsndnz_(j),this%ibsnd_(j,1:this%ibsndnz_(j)), &
                      this%ibsnddst_(j,1:this%ibsndnz_(j)),this%nbuff_,nt,j,v,this%sbbuff_)
    ENDDO
    DO j=1,this%ntsnd_  ! to top task:
      CALL gpc_packsf(nx,ny,nzl,this%itsndnz_(j),this%itsnd_(j,1:this%itsndnz_(j)), &
                      this%itsnddst_(j,1:this%itsndnz_(j)),this%nbuff_,nt,j,v,this%stbuff_)
    ENDDO

    CALL GTStart(this%hcomm_)
    CALL gpc_exch_planes(this,this%sbbuff_,this%stbuff_,this%rbbuff_,this%rtbuff_,nt)
    CALL GTAcc(this%hcomm_)

    DO j=1,this%nbrcv_
      CALL gpc_unpacksf(nx,ny,nex,ney,nez,ngp,this%ibrcvnz_(j),this%nbuff_,nt,j,this%rbbuff_,vext)
    ENDDO
    DO j=1,this%ntrcv_
      CALL gpc_unpacksf(nx,ny,nex,ney,nez,ngp,this%itrcvnz_(j),this%nbuff_,nt,j,this%rtbuff_,vext)
    ENDDO
  END SUBROUTINE GPartComm_SlabDataExchangeSF


! Posts the receives and sends of the ghost planes (device
! addresses of the buffers in offload builds) and waits for all
  SUBROUTINE gpc_exch_planes(this,sb,st,rb,rt,nt)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)               :: this
    INTEGER         ,INTENT(IN)                  :: nt
    REAL(KIND=GP),INTENT(IN)   ,TARGET :: sb(this%nbuff_,nt),st(this%nbuff_,nt)
    REAL(KIND=GP),INTENT(INOUT),TARGET :: rb(this%nbuff_,nt),rt(this%nbuff_,nt)
#if defined(GHOST_GPU)
    IF ( gdev_active ) THEN
!$omp target data use_device_addr(sb,st,rb,rt)
      CALL gpc_exch_planes_do(this,sb,st,rb,rt,nt)
!$omp end target data
    ELSE
      CALL gpc_exch_planes_do(this,sb,st,rb,rt,nt)
    ENDIF
#else
    CALL gpc_exch_planes_do(this,sb,st,rb,rt,nt)
#endif
  END SUBROUTINE gpc_exch_planes

  SUBROUTINE gpc_exch_planes_do(this,sb,st,rb,rt,nt)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)               :: this
    INTEGER         ,INTENT(IN)                  :: nt
    REAL(KIND=GP),INTENT(IN)    :: sb(this%nbuff_,nt),st(this%nbuff_,nt)
    REAL(KIND=GP),INTENT(INOUT) :: rb(this%nbuff_,nt),rt(this%nbuff_,nt)
    INTEGER                     :: j,buffsize
    buffsize = this%nd_(1)*this%nd_(2)*this%nzghost_ + this%nzghost_ + 2
    DO j=1,this%nbrcv_  ! from bottom task:
      CALL MPI_IRECV(rb(1,j),this%nbuff_,GC_REAL,this%ibrcvp_(j),1,this%comm_,this%ibrh_(j),this%ierr_)
    ENDDO
    DO j=1,this%ntrcv_  ! from top task:
      CALL MPI_IRECV(rt(1,j),this%nbuff_,GC_REAL,this%itrcvp_(j),2,this%comm_,this%itrh_(j),this%ierr_)
    ENDDO
    DO j=1,this%nbsnd_  ! to bottom task (arrives as its top data):
      CALL MPI_ISEND(sb(1,j),buffsize,GC_REAL,this%ibsndp_(j),2,this%comm_,this%ibsh_(j),this%ierr_)
    ENDDO
    DO j=1,this%ntsnd_  ! to top task (arrives as its bottom data):
      CALL MPI_ISEND(st(1,j),buffsize,GC_REAL,this%itsndp_(j),1,this%comm_,this%itsh_(j),this%ierr_)
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
  END SUBROUTINE gpc_exch_planes_do


! Interior of the extended field
  SUBROUTINE gpc_copy2ext(nx,ny,nzl,nex,ney,nez,ngp,ngz,v,vext)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: nx,ny,nzl,nex,ney,nez,ngp,ngz
    REAL(KIND=GP),INTENT(IN)    :: v(nx,ny,nzl)
    REAL(KIND=GP),INTENT(INOUT) :: vext(nex,ney,nez)
    INTEGER                     :: i,j,k
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private(i)
#endif
    DO k = 1,nzl
      DO j = 1,ny
        DO i = 1,nx
          vext(i+ngp,j+ngp,k+ngz) = v(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE gpc_copy2ext


! Ghost planes of a single task (periodic wrap of its own slab)
  SUBROUTINE gpc_localexch(nx,ny,nzl,nex,ney,nez,ngp,ngz,v,vext)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: nx,ny,nzl,nex,ney,nez,ngp,ngz
    REAL(KIND=GP),INTENT(IN)    :: v(nx,ny,nzl)
    REAL(KIND=GP),INTENT(INOUT) :: vext(nex,ney,nez)
    INTEGER                     :: i,j,k
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private(i)
#endif
    DO k = 1,ngz
      DO j = 1,ny
        DO i = 1,nx
          vext(i+ngp,j+ngp,nzl+ngz+k) = v(i,j,k)
          vext(i+ngp,j+ngp,k)         = v(i,j,nzl-ngz+k)
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE gpc_localexch


! Packs npl planes ksrc(1:npl) of v into column jc of buff, with
! the header (number of planes, destination plane indices kdst)
  SUBROUTINE gpc_packsf(nx,ny,nzl,npl,ksrc,kdst,nbuff,nt,jc,v,buff)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: nx,ny,nzl,npl,nbuff,nt,jc
    INTEGER      ,INTENT(IN)    :: ksrc(npl),kdst(npl)
    REAL(KIND=GP),INTENT(IN)    :: v(nx,ny,nzl)
    REAL(KIND=GP),INTENT(INOUT) :: buff(nbuff,nt)
    INTEGER                     :: i,j,m,nxy
    nxy = nx*ny
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
#endif
    DO m = 0,npl
      IF ( m.eq.0 ) THEN
        buff(1,jc) = REAL(npl,KIND=GP)
      ELSE
        buff(1+m,jc) = REAL(kdst(m),KIND=GP)
      ENDIF
    ENDDO
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private(i)
#endif
    DO m = 1,npl
      DO j = 1,ny
        DO i = 1,nx
          buff(1+npl+(m-1)*nxy+(j-1)*nx+i,jc) = v(i,j,ksrc(m))
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE gpc_packsf


! Unpacks the npl planes of column jc of buff into the extended
! field, at the plane indices carried by the header
  SUBROUTINE gpc_unpacksf(nx,ny,nex,ney,nez,ngp,npl,nbuff,nt,jc,buff,vext)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: nx,ny,nex,ney,nez,ngp,npl,nbuff,nt,jc
    REAL(KIND=GP),INTENT(IN)    :: buff(nbuff,nt)
    REAL(KIND=GP),INTENT(INOUT) :: vext(nex,ney,nez)
    INTEGER                     :: i,j,k,m,nxy
    nxy = nx*ny
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active) private(k)
#else
!$omp parallel do collapse(2) private(i,k)
#endif
    DO m = 1,npl
      DO j = 1,ny
        DO i = 1,nx
          k = INT(buff(1+m,jc))
          vext(i+ngp,j+ngp,k) = buff(1+npl+(m-1)*nxy+(j-1)*nx+i,jc)
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE gpc_unpacksf


!=================================================================
! Particle exchange between neighboring slabs
!=================================================================

!-----------------------------------------------------------------
!  METHOD     : IdentifyExchV
!  DESCRIPTION: Identifies the particles leaving the slab through
!               the bottom (pz < zmin) and through the top
!               (pz >= zmax): their indices are stored, in
!               ascending order, in ibot_(1:nbot_) and
!               itop_(1:ntop_). The counts are exchanged with the
!               neighbors and newnparts is the number of particles
!               the task will hold after the exchange. Positions
!               must not be periodized in z on entry.
!-----------------------------------------------------------------
  SUBROUTINE GPartComm_IdentifyExchV(this,id,pz,nparts,newnparts,zmin,zmax)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)           :: this
    INTEGER      ,INTENT   (IN)              :: nparts
    INTEGER      ,INTENT  (OUT)              :: newnparts
    INTEGER      ,INTENT(INOUT),DIMENSION(:) :: id
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:) :: pz
    REAL(KIND=GP),INTENT   (IN)              :: zmin,zmax
    INTEGER                                  :: ibrank,itrank,nrbot,nrtop,np

    newnparts = nparts
    IF ( this%nprocs_ .EQ. 1 ) RETURN ! nothing to do

    itrank = modulo(this%myrank_+1,this%nprocs_)
    ibrank = this%myrank_-1
    IF ( ibrank.LT.0 ) ibrank = this%nprocs_-1
    np = SIZE(pz)

    CALL gpc_flag_out(nparts,np,pz,zmin,zmax,1,this%iflag_)
    CALL gpsel_compact(1,nparts,this%iflag_,this%ibot_,this%nbot_,0)
    CALL gpc_flag_out(nparts,np,pz,zmin,zmax,0,this%iflag_)
    CALL gpsel_compact(1,nparts,this%iflag_,this%itop_,this%ntop_,0)

    CALL GTStart(this%hcomm_)
    CALL MPI_IRECV(nrbot,1,MPI_INTEGER,ibrank,1,this%comm_,this%ibrh_(1),this%ierr_)
    CALL MPI_IRECV(nrtop,1,MPI_INTEGER,itrank,2,this%comm_,this%itrh_(1),this%ierr_)
    CALL MPI_ISEND(this%nbot_,1,MPI_INTEGER,ibrank,2,this%comm_,this%ibsh_(1),this%ierr_)
    CALL MPI_ISEND(this%ntop_,1,MPI_INTEGER,itrank,1,this%comm_,this%itsh_(1),this%ierr_)
    CALL MPI_WAIT(this%ibrh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%itrh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%ibsh_(1),this%istatus_,this%ierr_)
    CALL MPI_WAIT(this%itsh_(1),this%istatus_,this%ierr_)
    CALL GTAcc(this%hcomm_)

    newnparts = nparts - (this%ntop_+this%nbot_) + (nrbot+nrtop)
  END SUBROUTINE GPartComm_IdentifyExchV


!-----------------------------------------------------------------
!  METHOD     : PartExchangeV
!  DESCRIPTION: Exchanges the particles identified by IdentifyExchV
!               with the neighbor tasks: the records (id and the
!               three arrays px,py,pz) leaving through the bottom
!               and the top are packed and sent, the departed
!               particles are removed and the received ones are
!               appended. With stg = GPEXCH_INIT or GPEXCH_UPDT
!               the id array and nparts are restored on exit, so
!               that the same lists can be used to exchange other
!               arrays of the same particles (positions of another
!               stage, velocities); GPEXCH_END or GPEXCH_UNIQ
!               commit the new list.
!-----------------------------------------------------------------
  SUBROUTINE GPartComm_PartExchangeV(this,id,px,py,pz,nparts,zmin,zmax,stg)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)           :: this
    INTEGER      ,INTENT(INOUT)              :: nparts
    INTEGER      ,INTENT(INOUT),DIMENSION(:) :: id
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:) :: px,py,pz
    REAL(KIND=GP),INTENT   (IN)              :: zmin,zmax
    INTEGER      ,INTENT   (IN), OPTIONAL    :: stg
    INTEGER                                  :: ibrank,itrank,ng,nrt,nrb,np,stage

    IF ( this%nprocs_ .EQ. 1 ) RETURN ! nothing to do

    itrank = modulo(this%myrank_+1,this%nprocs_)
    ibrank = this%myrank_-1
    IF ( ibrank.LT.0 ) ibrank = this%nprocs_-1
    np = SIZE(id)

    stage = GPEXCH_UNIQ
    IF (PRESENT(stg)) stage = stg
    ng = nparts
    IF ((stage.NE.GPEXCH_END).AND.(stage.NE.GPEXCH_UNIQ)) THEN
       CALL gpc_copy_i(ng,np,this%oldid_,id)
    ENDIF

    CALL gpc_pack(this%nbot_,np,this%ibot_,id,px,py,pz,this%sbid_,this%sbpr_)
    CALL gpc_pack(this%ntop_,np,this%itop_,id,px,py,pz,this%stid_,this%stpr_)
    CALL GTStart(this%hcomm_)
    CALL gpc_exch_parts(this,this%sbid_,this%sbpr_,this%stid_,this%stpr_, &
                        this%rbid_,this%rbpr_,this%rtid_,this%rtpr_,ibrank,itrank,nrb,nrt)
    CALL GTAcc(this%hcomm_)

    CALL GPartComm_ConcatV(this,id,px,py,pz,nparts)
    IF ( nparts+nrb+nrt .GT. np ) THEN
      WRITE(*,*) this%myrank_,' GPartComm_PartExchangeV: particle buffer too small: ', &
                 nparts+nrb+nrt,' > ',np
      STOP
    ENDIF
    CALL gpc_unpack(nrb,nparts,np,id,px,py,pz,this%rbid_,this%rbpr_)
    nparts = nparts + nrb
    CALL gpc_unpack(nrt,nparts,np,id,px,py,pz,this%rtid_,this%rtpr_)
    nparts = nparts + nrt

    IF ((stage.NE.GPEXCH_END).AND.(stage.NE.GPEXCH_UNIQ)) THEN
       CALL gpc_copy_i(ng,np,id,this%oldid_)
       nparts = ng
    ENDIF
  END SUBROUTINE GPartComm_PartExchangeV


! Sends the packed records to the two neighbors and receives
! theirs; nrb, nrt are the numbers of records received from the
! bottom and the top. Two messages per direction: ids, coordinates.
  SUBROUTINE gpc_exch_parts(this,sbid,sbpr,stid,stpr,rbid,rbpr,rtid,rtpr,ibrank,itrank,nrb,nrt)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)   :: this
    INTEGER      ,INTENT(IN)         :: ibrank,itrank
    INTEGER      ,INTENT(OUT)        :: nrb,nrt
    INTEGER      ,INTENT(IN)   ,TARGET :: sbid(this%maxparts_),stid(this%maxparts_)
    INTEGER      ,INTENT(INOUT),TARGET :: rbid(this%maxparts_),rtid(this%maxparts_)
    REAL(KIND=GP),INTENT(IN)   ,TARGET :: sbpr(3,this%maxparts_),stpr(3,this%maxparts_)
    REAL(KIND=GP),INTENT(INOUT),TARGET :: rbpr(3,this%maxparts_),rtpr(3,this%maxparts_)
#if defined(GHOST_GPU)
    IF ( gdev_active ) THEN
!$omp target data use_device_addr(sbid,sbpr,stid,stpr,rbid,rbpr,rtid,rtpr)
      CALL gpc_exch_parts_do(this,sbid,sbpr,stid,stpr,rbid,rbpr,rtid,rtpr,ibrank,itrank,nrb,nrt)
!$omp end target data
    ELSE
      CALL gpc_exch_parts_do(this,sbid,sbpr,stid,stpr,rbid,rbpr,rtid,rtpr,ibrank,itrank,nrb,nrt)
    ENDIF
#else
    CALL gpc_exch_parts_do(this,sbid,sbpr,stid,stpr,rbid,rbpr,rtid,rtpr,ibrank,itrank,nrb,nrt)
#endif
  END SUBROUTINE gpc_exch_parts

  SUBROUTINE gpc_exch_parts_do(this,sbid,sbpr,stid,stpr,rbid,rbpr,rtid,rtpr,ibrank,itrank,nrb,nrt)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)   :: this
    INTEGER      ,INTENT(IN)         :: ibrank,itrank
    INTEGER      ,INTENT(OUT)        :: nrb,nrt
    INTEGER      ,INTENT(IN)         :: sbid(*),stid(*)
    INTEGER      ,INTENT(INOUT)      :: rbid(*),rtid(*)
    REAL(KIND=GP),INTENT(IN)         :: sbpr(*),stpr(*)
    REAL(KIND=GP),INTENT(INOUT)      :: rbpr(*),rtpr(*)
    TYPE(MPI_Request)                :: req(4)
    TYPE(MPI_Status)                 :: st

    CALL MPI_ISEND(sbid,  this%nbot_,MPI_INTEGER,ibrank,1,this%comm_,req(1),this%ierr_)
    CALL MPI_ISEND(sbpr,3*this%nbot_,GC_REAL    ,ibrank,2,this%comm_,req(2),this%ierr_)
    CALL MPI_ISEND(stid,  this%ntop_,MPI_INTEGER,itrank,3,this%comm_,req(3),this%ierr_)
    CALL MPI_ISEND(stpr,3*this%ntop_,GC_REAL    ,itrank,4,this%comm_,req(4),this%ierr_)
    ! From the bottom neighbor come the records it sends to its top
    CALL MPI_RECV(rbid,  this%maxparts_,MPI_INTEGER,ibrank,3,this%comm_,st,this%ierr_)
    CALL MPI_GET_COUNT(st,MPI_INTEGER,nrb,this%ierr_)
    CALL MPI_RECV(rbpr,3*this%maxparts_,GC_REAL    ,ibrank,4,this%comm_,st,this%ierr_)
    CALL MPI_RECV(rtid,  this%maxparts_,MPI_INTEGER,itrank,1,this%comm_,st,this%ierr_)
    CALL MPI_GET_COUNT(st,MPI_INTEGER,nrt,this%ierr_)
    CALL MPI_RECV(rtpr,3*this%maxparts_,GC_REAL    ,itrank,2,this%comm_,st,this%ierr_)
    CALL MPI_WAITALL(4,req,MPI_STATUSES_IGNORE,this%ierr_)
  END SUBROUTINE gpc_exch_parts_do


!-----------------------------------------------------------------
!  METHOD     : ConcatV
!  DESCRIPTION: Removes the particles listed in ibot_ and itop_
!               and compacts the arrays: the holes left below the
!               new count are filled with the survivors found
!               above it, in ascending order on both sides, so the
!               result is deterministic (though not the original
!               order). nparts is updated.
!-----------------------------------------------------------------
  SUBROUTINE GPartComm_ConcatV(this,id,px,py,pz,nparts)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)           :: this
    INTEGER      ,INTENT(INOUT)              :: nparts
    INTEGER      ,INTENT(INOUT),DIMENSION(:) :: id
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:) :: px,py,pz
    INTEGER                                  :: ngood,nh,ns,np

    IF ((this%nbot_+this%ntop_).EQ.0) RETURN ! nothing to do
    np    = SIZE(id)
    ngood = nparts - (this%nbot_+this%ntop_)
    CALL gpc_mark(this%nbot_,np,this%ibot_,id)
    CALL gpc_mark(this%ntop_,np,this%itop_,id)
    CALL gpc_flag_null(1,ngood,np,id,.TRUE.,this%iflag_)
    CALL gpsel_compact(1,ngood,this%iflag_,this%ihole_,nh,0)
    CALL gpc_flag_null(ngood+1,nparts,np,id,.FALSE.,this%iflag_)
    CALL gpsel_compact(ngood+1,nparts,this%iflag_,this%isurv_,ns,0)
    IF ( nh .NE. ns ) THEN
      WRITE(*,*) this%myrank_,' GPartComm_ConcatV: inconsistent compaction: ',nh,ns
      STOP
    ENDIF
    CALL gpc_move(nh,np,this%ihole_,this%isurv_,id,px,py,pz)
    nparts = ngood
  END SUBROUTINE GPartComm_ConcatV


!-----------------------------------------------------------------
! Particle kernels (explicit-shape arrays; np is the size of the
! particle arrays, n the number of entries to process)
!-----------------------------------------------------------------

! flag(j) = 1 for the particles below zmin (ibelow=1) or at/above
! zmax (ibelow=0), 0 otherwise
  SUBROUTINE gpc_flag_out(n,np,pz,zmin,zmax,ibelow,flag)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n,np,ibelow
    REAL(KIND=GP),INTENT(IN)    :: pz(np)
    REAL(KIND=GP),INTENT(IN)    :: zmin,zmax
    INTEGER      ,INTENT(INOUT) :: flag(np)
    INTEGER                     :: j
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
#else
!$omp parallel do
#endif
    DO j = 1,n
      IF ( ibelow.eq.1 ) THEN
        IF ( pz(j).LT.zmin ) THEN
          flag(j) = 1
        ELSE
          flag(j) = 0
        ENDIF
      ELSE
        IF ( pz(j).GE.zmax ) THEN
          flag(j) = 1
        ELSE
          flag(j) = 0
        ENDIF
      ENDIF
    ENDDO
  END SUBROUTINE gpc_flag_out


! Gathers the records idx(1:n) into the send buffers
  SUBROUTINE gpc_pack(n,np,idx,id,px,py,pz,bid,bpr)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n,np
    INTEGER      ,INTENT(IN)    :: idx(np),id(np)
    REAL(KIND=GP),INTENT(IN)    :: px(np),py(np),pz(np)
    INTEGER      ,INTENT(INOUT) :: bid(np)
    REAL(KIND=GP),INTENT(INOUT) :: bpr(3,np)
    INTEGER                     :: j,i
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active) private(i)
#else
!$omp parallel do private(i)
#endif
    DO j = 1,n
      i = idx(j)
      bid(j)   = id(i)
      bpr(1,j) = px(i)
      bpr(2,j) = py(i)
      bpr(3,j) = pz(i)
    ENDDO
  END SUBROUTINE gpc_pack


! Appends the n received records after entry ioff
  SUBROUTINE gpc_unpack(n,ioff,np,id,px,py,pz,bid,bpr)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n,ioff,np
    INTEGER      ,INTENT(INOUT) :: id(np)
    REAL(KIND=GP),INTENT(INOUT) :: px(np),py(np),pz(np)
    INTEGER      ,INTENT(IN)    :: bid(np)
    REAL(KIND=GP),INTENT(IN)    :: bpr(3,np)
    INTEGER                     :: j
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
#else
!$omp parallel do
#endif
    DO j = 1,n
      id(ioff+j) = bid(j)
      px(ioff+j) = bpr(1,j)
      py(ioff+j) = bpr(2,j)
      pz(ioff+j) = bpr(3,j)
    ENDDO
  END SUBROUTINE gpc_unpack


! Marks the particles idx(1:n) as departed
  SUBROUTINE gpc_mark(n,np,idx,id)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n,np
    INTEGER      ,INTENT(IN)    :: idx(np)
    INTEGER      ,INTENT(INOUT) :: id(np)
    INTEGER                     :: j
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
#else
!$omp parallel do
#endif
    DO j = 1,n
      id(idx(j)) = GPNULL
    ENDDO
  END SUBROUTINE gpc_mark


! flag(j) = 1 where (id(j) == GPNULL) .eqv. wantnull, for j in [n1,n2]
  SUBROUTINE gpc_flag_null(n1,n2,np,id,wantnull,flag)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n1,n2,np
    INTEGER      ,INTENT(IN)    :: id(np)
    LOGICAL      ,INTENT(IN)    :: wantnull
    INTEGER      ,INTENT(INOUT) :: flag(np)
    INTEGER                     :: j
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
#else
!$omp parallel do
#endif
    DO j = n1,n2
      IF ( (id(j).EQ.GPNULL) .EQV. wantnull ) THEN
        flag(j) = 1
      ELSE
        flag(j) = 0
      ENDIF
    ENDDO
  END SUBROUTINE gpc_flag_null


! Fills the holes ihole(1:m) with the survivors isurv(1:m)
  SUBROUTINE gpc_move(m,np,ihole,isurv,id,px,py,pz)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: m,np
    INTEGER      ,INTENT(IN)    :: ihole(np),isurv(np)
    INTEGER      ,INTENT(INOUT) :: id(np)
    REAL(KIND=GP),INTENT(INOUT) :: px(np),py(np),pz(np)
    INTEGER                     :: k,ih,is
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active) private(ih,is)
#else
!$omp parallel do private(ih,is)
#endif
    DO k = 1,m
      ih = ihole(k)
      is = isurv(k)
      id(ih) = id(is)
      px(ih) = px(is)
      py(ih) = py(is)
      pz(ih) = pz(is)
    ENDDO
  END SUBROUTINE gpc_move


! dst(1:n) = src(1:n)
  SUBROUTINE gpc_copy_i(n,np,dst,src)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n,np
    INTEGER      ,INTENT(INOUT) :: dst(np)
    INTEGER      ,INTENT(IN)    :: src(np)
    INTEGER                     :: j
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
#else
!$omp parallel do
#endif
    DO j = 1,n
      dst(j) = src(j)
    ENDDO
  END SUBROUTINE gpc_copy_i


! idx(j) = j
  SUBROUTINE gpc_iota(n,np,idx)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n,np
    INTEGER      ,INTENT(INOUT) :: idx(np)
    INTEGER                     :: j
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
#else
!$omp parallel do
#endif
    DO j = 1,n
      idx(j) = j
    ENDDO
  END SUBROUTINE gpc_iota


! a(1:3,1:n) = 0
  SUBROUTINE gpc_zero3(n,a)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n
    REAL(KIND=GP),INTENT(INOUT) :: a(3,n)
    INTEGER                     :: j
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
#else
!$omp parallel do
#endif
    DO j = 1,n
      a(1,j) = 0.0_GP
      a(2,j) = 0.0_GP
      a(3,j) = 0.0_GP
    ENDDO
  END SUBROUTINE gpc_zero3


! Scatters the nl local records into the global array by id
  SUBROUTINE gpc_scatter3(nl,np,id,lx,ly,lz,ng,g)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: nl,np,ng
    INTEGER      ,INTENT(IN)    :: id(np)
    REAL(KIND=GP),INTENT(IN)    :: lx(np),ly(np),lz(np)
    REAL(KIND=GP),INTENT(INOUT) :: g(3,ng)
    INTEGER                     :: j,i
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active) private(i)
#else
!$omp parallel do private(i)
#endif
    DO j = 1,nl
      i = id(j) + 1
      g(1,i) = lx(j)
      g(2,i) = ly(j)
      g(3,i) = lz(j)
    ENDDO
  END SUBROUTINE gpc_scatter3


!=================================================================
! Global particle database (VDB exchange)
!=================================================================

!-----------------------------------------------------------------
!  METHOD     : VDBSynch
!  DESCRIPTION: Builds the global database gvdb(3,ngvdb) of the
!               positions of all the particles, indexed by id,
!               from the local records of every task (scatter by
!               id into ptmp, then a sum over tasks).
!-----------------------------------------------------------------
  SUBROUTINE GPartComm_VDBSynch(this,gvdb,ngvdb,id,lx,ly,lz,nl,ptmp)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)                       :: this
    INTEGER      ,INTENT   (IN),DIMENSION(:)             :: id
    INTEGER      ,INTENT   (IN)                          :: nl
    INTEGER      ,INTENT   (IN)                          :: ngvdb
    REAL(KIND=GP),INTENT   (IN),DIMENSION(:)             :: lx,ly,lz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(3,ngvdb),TARGET:: gvdb,ptmp

    CALL gpc_zero3(ngvdb,ptmp)
    CALL gpc_scatter3(nl,SIZE(id),id,lx,ly,lz,ngvdb,ptmp)
    CALL GTStart(this%hcomm_)
#if defined(GHOST_GPU)
    IF ( gdev_active ) THEN
!$omp target data use_device_addr(ptmp,gvdb)
      CALL MPI_ALLREDUCE(ptmp,gvdb,3*ngvdb,GC_REAL,MPI_SUM,this%comm_,this%ierr_)
!$omp end target data
    ELSE
      CALL MPI_ALLREDUCE(ptmp,gvdb,3*ngvdb,GC_REAL,MPI_SUM,this%comm_,this%ierr_)
    ENDIF
#else
    CALL MPI_ALLREDUCE(ptmp,gvdb,3*ngvdb,GC_REAL,MPI_SUM,this%comm_,this%ierr_)
#endif
    CALL GTAcc(this%hcomm_)
  END SUBROUTINE GPartComm_VDBSynch


!-----------------------------------------------------------------
!  METHOD     : VDBSynch_t0
!  DESCRIPTION: Gathers on task 0 the records of all the tasks
!               into gvdb, indexed by id (I/O; host copies)
!-----------------------------------------------------------------
  SUBROUTINE GPartComm_VDBSynch_t0(this,gvdb,ngvdb,id,lx,ly,lz,nl)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)                  :: this
    INTEGER      ,INTENT   (IN),DIMENSION(:)        :: id
    INTEGER      ,INTENT   (IN)                     :: nl
    INTEGER      ,INTENT   (IN)                     :: ngvdb
    REAL(KIND=GP),INTENT   (IN),DIMENSION(:)        :: lx,ly,lz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(3,ngvdb)  :: gvdb
    INTEGER                                         :: t,i,j,m,np
    TYPE(MPI_Request)                               :: req(2)

    IF ( gdev_active ) STOP 'GPartComm_VDBSynch_t0: must be called with gdev_active unset'
    np = SIZE(id)
    IF (this%myrank_.NE.0) THEN
      CALL gpc_iota(nl,np,this%itop_)
      CALL gpc_pack(nl,np,this%itop_,id,lx,ly,lz,this%stid_,this%stpr_)
      CALL GTStart(this%hcomm_)
      CALL MPI_ISEND(this%stid_,  nl,MPI_INTEGER,0,this%myrank_             ,this%comm_,req(1),this%ierr_)
      CALL MPI_ISEND(this%stpr_,3*nl,GC_REAL    ,0,this%myrank_+this%nprocs_,this%comm_,req(2),this%ierr_)
      CALL MPI_WAITALL(2,req,MPI_STATUSES_IGNORE,this%ierr_)
      CALL GTAcc(this%hcomm_)
    ELSE
      DO j = 1,nl
        i = id(j) + 1
        gvdb(1,i) = lx(j)
        gvdb(2,i) = ly(j)
        gvdb(3,i) = lz(j)
      END DO
      DO t = 1,this%nprocs_-1
        CALL GTStart(this%hcomm_)
        CALL MPI_RECV(this%rbid_,  this%maxparts_,MPI_INTEGER,t,t             ,this%comm_,this%istatus_,this%ierr_)
        CALL MPI_GET_COUNT(this%istatus_,MPI_INTEGER,m,this%ierr_)
        CALL MPI_RECV(this%rbpr_,3*this%maxparts_,GC_REAL    ,t,t+this%nprocs_,this%comm_,this%istatus_,this%ierr_)
        CALL GTAcc(this%hcomm_)
        DO j = 1,m
          i = this%rbid_(j) + 1
          gvdb(1,i) = this%rbpr_(1,j)
          gvdb(2,i) = this%rbpr_(2,j)
          gvdb(3,i) = this%rbpr_(3,j)
        END DO
      END DO
    END IF
  END SUBROUTINE GPartComm_VDBSynch_t0


!-----------------------------------------------------------------
!  METHOD     : LagSynch_t0
!  DESCRIPTION: Gathers on task 0 a scalar Lagrangian quantity ls
!               of all the tasks into gs, indexed by id (I/O)
!-----------------------------------------------------------------
  SUBROUTINE GPartComm_LagSynch_t0(this,gs,ngs,id,ls,nl)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)                  :: this
    INTEGER      ,INTENT   (IN),DIMENSION(:)        :: id
    INTEGER      ,INTENT   (IN)                     :: nl
    INTEGER      ,INTENT   (IN)                     :: ngs
    REAL(KIND=GP),INTENT   (IN),DIMENSION(:)        :: ls
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:)        :: gs
    INTEGER                                         :: i,j,t,m
    TYPE(MPI_Request)                               :: req(2)

    IF ( gdev_active ) STOP 'GPartComm_LagSynch_t0: must be called with gdev_active unset'
    IF (this%myrank_.NE.0) THEN
      DO i = 1,nl
        this%stid_(i)   = id(i)
        this%stpr_(1,i) = ls(i)
      END DO 
      CALL GTStart(this%hcomm_)
      CALL MPI_ISEND(this%stid_,  nl,MPI_INTEGER,0,this%myrank_             ,this%comm_,req(1),this%ierr_)
      CALL MPI_ISEND(this%stpr_,3*nl,GC_REAL    ,0,this%myrank_+this%nprocs_,this%comm_,req(2),this%ierr_)
      CALL MPI_WAITALL(2,req,MPI_STATUSES_IGNORE,this%ierr_)
      CALL GTAcc(this%hcomm_)
    ELSE
      DO j = 1,nl
        i = id(j) + 1
        gs(i) = ls(j)
      END DO
      DO t = 1,this%nprocs_-1
        CALL GTStart(this%hcomm_)
        CALL MPI_RECV(this%rbid_,  this%maxparts_,MPI_INTEGER,t,t             ,this%comm_,this%istatus_,this%ierr_)
        CALL MPI_GET_COUNT(this%istatus_,MPI_INTEGER,m,this%ierr_)
        CALL MPI_RECV(this%rbpr_,3*this%maxparts_,GC_REAL    ,t,t+this%nprocs_,this%comm_,this%istatus_,this%ierr_)
        CALL GTAcc(this%hcomm_)
        DO j = 1,m
          i     = this%rbid_(j) + 1
          gs(i) = this%rbpr_(1,j)
        END DO
      END DO
    END IF
  END SUBROUTINE GPartComm_LagSynch_t0


!=================================================================
! Initial distribution of the particles read by task 0
!=================================================================

!-----------------------------------------------------------------
!  METHOD     : IdentifyTaskV
!  DESCRIPTION: Task 0 finds the task owning each particle from
!               its z coordinate (task(i)) and tells every task how
!               many particles it will receive (host, at init)
!-----------------------------------------------------------------
  SUBROUTINE GPartComm_IdentifyTaskV(this,id,pz,nparts,task)
    USE grid
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)           :: this
    INTEGER      ,INTENT(INOUT)              :: nparts
    INTEGER                                  :: i,j,tsta,tend
    INTEGER      ,INTENT(INOUT),DIMENSION(:) :: id,task
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:) :: pz
    REAL(KIND=GP)                            :: zmaxs(this%nprocs_)
    INTEGER                                  :: nps(this%nprocs_)
    TYPE(MPI_Request)                        :: req(this%nprocs_)

    IF ( this%nprocs_ .EQ. 1 ) RETURN ! nothing to do
    IF (this%myrank_.EQ.0) THEN
      DO j = 1,this%nprocs_
        nps(j) = 0
        CALL range(1,nz,this%nprocs_,j-1,tsta,tend)
        zmaxs(j) = REAL(tend,KIND=GP)
      END DO
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
      DO j = 2,this%nprocs_
        CALL MPI_ISEND(nps(j),1,MPI_INTEGER,j-1,j-1,this%comm_,req(j),this%ierr_)
      END DO
      DO j = 2,this%nprocs_
        CALL MPI_WAIT(req(j),this%istatus_,this%ierr_)
      END DO
    ELSE
      CALL MPI_RECV(nparts,1,MPI_INTEGER,0,this%myrank_,this%comm_,this%istatus_,this%ierr_)
    END IF
  END SUBROUTINE GPartComm_IdentifyTaskV


!-----------------------------------------------------------------
!  METHOD     : PartScatterV
!  DESCRIPTION: Task 0 sends each task the particles it owns
!               according to task(:) and keeps its own (host, at
!               init; the callers upload the states afterwards)
!-----------------------------------------------------------------
  SUBROUTINE GPartComm_PartScatterV(this,id,px,py,pz,nparts,task)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)           :: this
    INTEGER      ,INTENT(INOUT)              :: nparts
    INTEGER                                  :: i,j,nsend,nparts0,np
    INTEGER      ,INTENT(INOUT),DIMENSION(:) :: id
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:) :: px,py,pz
    INTEGER      ,INTENT   (IN),DIMENSION(:) :: task
    TYPE(MPI_Request)                        :: req(2)

    IF ( this%nprocs_ .EQ. 1 ) RETURN ! nothing to do
    IF ( gdev_active ) STOP 'GPartComm_PartScatterV: must be called with gdev_active unset'
    np = SIZE(id)
    IF (this%myrank_.EQ.0) THEN
      DO j = 2,this%nprocs_
        nsend = 0
        DO i = 1,this%maxparts_
          IF (task(i).EQ.(j-1)) THEN
            nsend = nsend + 1
            this%itop_(nsend) = i
          END IF
        END DO
        CALL gpc_pack(nsend,np,this%itop_,id,px,py,pz,this%stid_,this%stpr_)
        CALL MPI_ISEND(this%stid_,  nsend,MPI_INTEGER,j-1,j-1             ,this%comm_,req(1),this%ierr_)
        CALL MPI_ISEND(this%stpr_,3*nsend,GC_REAL    ,j-1,j-1+this%nprocs_,this%comm_,req(2),this%ierr_)
        CALL MPI_WAITALL(2,req,MPI_STATUSES_IGNORE,this%ierr_)
      END DO
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
      CALL MPI_RECV(this%rtid_,  this%maxparts_,MPI_INTEGER,0,this%myrank_             ,this%comm_,this%istatus_,this%ierr_)
      CALL MPI_GET_COUNT(this%istatus_,MPI_INTEGER,nparts,this%ierr_)
      CALL MPI_RECV(this%rtpr_,3*this%maxparts_,GC_REAL    ,0,this%myrank_+this%nprocs_,this%comm_,this%istatus_,this%ierr_)
      CALL gpc_unpack(nparts,0,np,id,px,py,pz,this%rtid_,this%rtpr_)
    END IF
  END SUBROUTINE GPartComm_PartScatterV


!-----------------------------------------------------------------
!  METHOD     : ResizeArrays
!  DESCRIPTION: Resizes the particle-sized arrays to newmparts
!               (only grows them if onlyinc). The exchange lists
!               keep their contents.
!-----------------------------------------------------------------
  SUBROUTINE GPartComm_ResizeArrays(this,newmparts,onlyinc)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)  :: this
    INTEGER         ,INTENT(IN)     :: newmparts
    LOGICAL         ,INTENT(IN)     :: onlyinc
    INTEGER                         :: n

    n = SIZE(this%ibot_)
    IF ((n.lt.newmparts).OR.((n.gt.newmparts).AND..NOT.onlyinc)) THEN
      CALL gresize(this%ibot_ ,newmparts,.true. )
      CALL gresize(this%itop_ ,newmparts,.true. )
      CALL gresize(this%oldid_,newmparts,.false.)
      CALL gresize(this%iflag_,newmparts,.false.)
      CALL gresize(this%ihole_,newmparts,.false.)
      CALL gresize(this%isurv_,newmparts,.false.)
      CALL gresize(this%sbid_ ,newmparts,.false.)
      CALL gresize(this%stid_ ,newmparts,.false.)
      CALL gresize(this%rbid_ ,newmparts,.false.)
      CALL gresize(this%rtid_ ,newmparts,.false.)
      CALL gresize(this%sbpr_ ,3,newmparts,.false.)
      CALL gresize(this%stpr_ ,3,newmparts,.false.)
      CALL gresize(this%rbpr_ ,3,newmparts,.false.)
      CALL gresize(this%rtpr_ ,3,newmparts,.false.)
    ENDIF
    this%maxparts_ = newmparts
  END SUBROUTINE GPartComm_ResizeArrays


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
    ! Field buffers, particle buffers and index lists have device copies
    CALL galloc(this%sbbuff_,this%nbuff_,nt)
    CALL galloc(this%stbuff_,this%nbuff_,nt)
    CALL galloc(this%rbbuff_,this%nbuff_,nt)
    CALL galloc(this%rtbuff_,this%nbuff_,nt)
    CALL GPartComm_AllocParts(this,this%maxparts_)
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
            this%itsnd_   (this%ntsnd_,k) = kend-ksta+2-k                  ! local z-index to be sent to top
            this%itsnddst_(this%ntsnd_,k) = this%nzghost_-kf-k+1           ! local destination z-index
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
            this%ibsnd_   (this%nbsnd_,k) = k    ! local z-index to be sent to bottom
            this%ibsnddst_(this%nbsnd_,k) = this%nzghost_+nzf(jf)+kf+k     ! local destination z-index
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

    DEALLOCATE(jfwd,kfend,kfsta,nzF)

  END SUBROUTINE GPartComm_Init

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

!=================================================================
! Transposes of a real field between the slab layout (nx,ny,kl),
! kl = kend-ksta+1 local planes, and the z-complete layout
! (nz,ny,il), il = local x range, used by the spline solve in z.
! The blocks destined to each task are packed into a contiguous
! device buffer, exchanged with one message per task pair, and
! unpacked on the device (GPU-aware MPI; the buffers keep their
! device copies between calls). The block sent to task t in the
! forward transpose has the same layout as the block received
! from t in the inverse one, so the offsets and counts are shared.
!=================================================================

! Tables of the exchange: x and z ranges of every task, offsets
! and counts (in reals) of the blocks; buffers of the size of the
! larger of the two layouts
  SUBROUTINE GPartComm_InitTrans(this)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT) :: this
    INTEGER                        :: t,nx,ny,nz,kl,il,ntot,ttot,nb

    nx = this%nd_(1); ny = this%nd_(2); nz = this%nd_(3)
    ALLOCATE(this%txsta_(0:this%nprocs_-1),this%txend_(0:this%nprocs_-1))
    ALLOCATE(this%tzsta_(0:this%nprocs_-1),this%tzend_(0:this%nprocs_-1))
    ALLOCATE(this%tso_(0:this%nprocs_-1),this%tsc_(0:this%nprocs_-1))
    ALLOCATE(this%tro_(0:this%nprocs_-1),this%trc_(0:this%nprocs_-1))
    DO t = 0,this%nprocs_-1
      CALL range(1,nx,this%nprocs_,t,this%txsta_(t),this%txend_(t))
      CALL range(1,nz,this%nprocs_,t,this%tzsta_(t),this%tzend_(t))
    ENDDO
    kl = this%kend_-this%ksta_+1
    il = this%txend_(this%myrank_)-this%txsta_(this%myrank_)+1
    ntot = 0; ttot = 0
    DO t = 0,this%nprocs_-1
      this%tso_(t) = ntot                    ! block (x range of t, ny, my kl planes)
      this%tsc_(t) = (this%txend_(t)-this%txsta_(t)+1)*ny*kl
      ntot = ntot + this%tsc_(t)
      this%tro_(t) = ttot                    ! block (my x range, ny, z range of t)
      this%trc_(t) = il*ny*(this%tzend_(t)-this%tzsta_(t)+1)
      ttot = ttot + this%trc_(t)
    ENDDO
    nb = MAX(ntot,ttot)
    CALL galloc(this%gtsbuf_,nb)
    CALL galloc(this%gtrbuf_,nb)
    this%btransinit_ = .TRUE.
  END SUBROUTINE GPartComm_InitTrans


!-----------------------------------------------------------------
!  METHOD     : GTranspose
!  DESCRIPTION: ofield(nz,ny,il) = transpose of ifield(nx,ny,kl)
!-----------------------------------------------------------------
  SUBROUTINE GPartComm_Transpose(this,ofield,ifield)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)           :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*) :: ofield
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*) :: ifield
    INTEGER                                  :: t,nx,ny,nz,kl,il,ista,ksta

    IF ( .NOT.this%btransinit_ ) CALL GPartComm_InitTrans(this)
    nx = this%nd_(1); ny = this%nd_(2); nz = this%nd_(3)
    kl = this%kend_-this%ksta_+1
    ista = this%txsta_(this%myrank_); il = this%txend_(this%myrank_)-ista+1
    ksta = this%ksta_
    DO t = 0,this%nprocs_-1
      CALL gpc_tpack_fwd(nx,ny,kl,this%txsta_(t),this%txend_(t),this%tso_(t),SIZE(this%gtsbuf_),ifield,this%gtsbuf_)
    ENDDO
    CALL GTStart(this%hcomm_)
    CALL gpc_texch(this,SIZE(this%gtsbuf_),this%gtsbuf_,this%gtrbuf_,this%tso_,this%tsc_,this%tro_,this%trc_)
    CALL GTAcc(this%hcomm_)
    DO t = 0,this%nprocs_-1
      CALL gpc_tunpack_fwd(nz,ny,il,this%tzsta_(t),this%tzend_(t),this%tro_(t),SIZE(this%gtrbuf_),this%gtrbuf_,ofield)
    ENDDO
  END SUBROUTINE GPartComm_Transpose


!-----------------------------------------------------------------
!  METHOD     : GITranspose
!  DESCRIPTION: ofield(nx,ny,kl) = inverse transpose of ifield(nz,ny,il)
!-----------------------------------------------------------------
  SUBROUTINE GPartComm_ITranspose(this,ofield,ifield)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)           :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(*) :: ofield
    REAL(KIND=GP),INTENT   (IN),DIMENSION(*) :: ifield
    INTEGER                                  :: t,nx,ny,nz,kl,il,ista

    IF ( .NOT.this%btransinit_ ) CALL GPartComm_InitTrans(this)
    nx = this%nd_(1); ny = this%nd_(2); nz = this%nd_(3)
    kl = this%kend_-this%ksta_+1
    ista = this%txsta_(this%myrank_); il = this%txend_(this%myrank_)-ista+1
    DO t = 0,this%nprocs_-1
      CALL gpc_tpack_inv(nz,ny,il,this%tzsta_(t),this%tzend_(t),this%tro_(t),SIZE(this%gtsbuf_),ifield,this%gtsbuf_)
    ENDDO
    CALL GTStart(this%hcomm_)
    CALL gpc_texch(this,SIZE(this%gtsbuf_),this%gtsbuf_,this%gtrbuf_,this%tro_,this%trc_,this%tso_,this%tsc_)
    CALL GTAcc(this%hcomm_)
    DO t = 0,this%nprocs_-1
      CALL gpc_tunpack_inv(nx,ny,kl,this%txsta_(t),this%txend_(t),this%tso_(t),SIZE(this%gtrbuf_),this%gtrbuf_,ofield)
    ENDDO
  END SUBROUTINE GPartComm_ITranspose


! Contiguous exchange of the blocks (device addresses in offload
! builds): the block for task t starts at so(t)+1 in sb and the
! block from t lands at ro(t)+1 in rb; the own block is copied
  SUBROUTINE gpc_texch(this,nb,sb,rb,so,sc,ro,rc)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)   :: this
    INTEGER      ,INTENT(IN)           :: nb
    ! Explicit-shape buffers: with assumed-shape dummies the compiler
    ! may copy them through the host when they are passed on inside
    ! the use_device_addr region, reading device memory from the CPU
    REAL(KIND=GP),INTENT(IN)   ,TARGET :: sb(nb)
    REAL(KIND=GP),INTENT(INOUT),TARGET :: rb(nb)
    INTEGER      ,INTENT(IN)           :: so(0:this%nprocs_-1),sc(0:this%nprocs_-1)
    INTEGER      ,INTENT(IN)           :: ro(0:this%nprocs_-1),rc(0:this%nprocs_-1)
    CALL gpc_copy_seg(sc(this%myrank_),nb,so(this%myrank_),ro(this%myrank_),sb,rb)
    IF ( this%nprocs_ .EQ. 1 ) RETURN
#if defined(GHOST_GPU)
    IF ( gdev_active ) THEN
!$omp target data use_device_addr(sb,rb)
      CALL gpc_texch_do(this,sb,rb,so,sc,ro,rc)
!$omp end target data
    ELSE
      CALL gpc_texch_do(this,sb,rb,so,sc,ro,rc)
    ENDIF
#else
    CALL gpc_texch_do(this,sb,rb,so,sc,ro,rc)
#endif
  END SUBROUTINE gpc_texch

  SUBROUTINE gpc_texch_do(this,sb,rb,so,sc,ro,rc)
    IMPLICIT NONE
    CLASS(GPartComm),INTENT(INOUT)   :: this
    REAL(KIND=GP),INTENT(IN)         :: sb(*)
    REAL(KIND=GP),INTENT(INOUT)      :: rb(*)
    INTEGER      ,INTENT(IN)         :: so(0:this%nprocs_-1),sc(0:this%nprocs_-1)
    INTEGER      ,INTENT(IN)         :: ro(0:this%nprocs_-1),rc(0:this%nprocs_-1)
    TYPE(MPI_Request)                :: req(2*this%nprocs_)
    INTEGER                          :: irank,isendto,igetfrom,nreq
    nreq = 0
    DO irank = 1,this%nprocs_-1
      igetfrom = modulo(this%myrank_-irank+this%nprocs_,this%nprocs_)
      nreq = nreq+1
      CALL MPI_IRECV(rb(ro(igetfrom)+1),rc(igetfrom),GC_REAL,igetfrom,1,this%comm_,req(nreq),this%ierr_)
    ENDDO
    DO irank = 1,this%nprocs_-1
      isendto = modulo(this%myrank_+irank,this%nprocs_)
      nreq = nreq+1
      CALL MPI_ISEND(sb(so(isendto)+1),sc(isendto),GC_REAL,isendto,1,this%comm_,req(nreq),this%ierr_)
    ENDDO
    CALL MPI_WAITALL(nreq,req(1:nreq),MPI_STATUSES_IGNORE,this%ierr_)
  END SUBROUTINE gpc_texch_do


! rb(roff+1:roff+n) = sb(soff+1:soff+n)
  SUBROUTINE gpc_copy_seg(n,nb,soff,roff,sb,rb)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: n,nb,soff,roff
    REAL(KIND=GP),INTENT(IN)    :: sb(nb)
    REAL(KIND=GP),INTENT(INOUT) :: rb(nb)
    INTEGER                     :: i
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
#else
!$omp parallel do
#endif
    DO i = 1,n
      rb(roff+i) = sb(soff+i)
    ENDDO
  END SUBROUTINE gpc_copy_seg


! Forward pack: block f(i1:i2,1:ny,1:kl) of the slab, i fastest
  SUBROUTINE gpc_tpack_fwd(nx,ny,kl,i1,i2,off,nb,f,b)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: nx,ny,kl,i1,i2,off,nb
    REAL(KIND=GP),INTENT(IN)    :: f(nx,ny,kl)
    REAL(KIND=GP),INTENT(INOUT) :: b(nb)
    INTEGER                     :: i,j,k,ni
    ni = i2-i1+1
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private(i)
#endif
    DO k = 1,kl
      DO j = 1,ny
        DO i = i1,i2
          b(off+(i-i1+1)+(j-1)*ni+(k-1)*ni*ny) = f(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE gpc_tpack_fwd


! Forward unpack: block (1:il,1:ny,k1:k2) received from a task
! into the z-complete layout o(nz,ny,il)
  SUBROUTINE gpc_tunpack_fwd(nz,ny,il,k1,k2,off,nb,b,o)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: nz,ny,il,k1,k2,off,nb
    REAL(KIND=GP),INTENT(IN)    :: b(nb)
    REAL(KIND=GP),INTENT(INOUT) :: o(nz,ny,il)
    INTEGER                     :: i,j,k
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private(k)
#endif
    DO i = 1,il
      DO j = 1,ny
        DO k = k1,k2
          o(k,j,i) = b(off+i+(j-1)*il+(k-k1)*il*ny)
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE gpc_tunpack_fwd


! Inverse pack: block (k1:k2,1:ny,1:il) of the z-complete layout,
! stored i fastest (the layout the receiver unpacks)
  SUBROUTINE gpc_tpack_inv(nz,ny,il,k1,k2,off,nb,f,b)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: nz,ny,il,k1,k2,off,nb
    REAL(KIND=GP),INTENT(IN)    :: f(nz,ny,il)
    REAL(KIND=GP),INTENT(INOUT) :: b(nb)
    INTEGER                     :: i,j,k
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private(k)
#endif
    DO i = 1,il
      DO j = 1,ny
        DO k = k1,k2
          b(off+i+(j-1)*il+(k-k1)*il*ny) = f(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE gpc_tpack_inv


! Inverse unpack: block (i1:i2,1:ny,1:kl) into the slab o(nx,ny,kl)
  SUBROUTINE gpc_tunpack_inv(nx,ny,kl,i1,i2,off,nb,b,o)
    IMPLICIT NONE
    INTEGER      ,INTENT(IN)    :: nx,ny,kl,i1,i2,off,nb
    REAL(KIND=GP),INTENT(IN)    :: b(nb)
    REAL(KIND=GP),INTENT(INOUT) :: o(nx,ny,kl)
    INTEGER                     :: i,j,k,ni
    ni = i2-i1+1
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private(i)
#endif
    DO k = 1,kl
      DO j = 1,ny
        DO i = i1,i2
          o(i,j,k) = b(off+(i-i1+1)+(j-1)*ni+(k-1)*ni*ny)
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE gpc_tunpack_inv


END MODULE class_GPartComm
