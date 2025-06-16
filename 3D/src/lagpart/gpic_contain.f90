!=================================================================
! GPIC SUBROUTINES
!=================================================================

  SUBROUTINE GPIC_ctor(this,comm,ppc,inittype,intorder,iexchtyp, &
                        iouttyp,bcoll,csize,nstrip,intacc,wrtunit)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Explicit constructor for particle-in-cell. Should be called
!  after calling GPart_ctor.
!
!  ARGUMENTS:
!    this    : 'this' class instance
!    comm    : MPI communicator
!    ppc     : number of particle per cell
!    inittype: GPINIT-typed quantity to give the type of particle
!              position initialization
!    intorder: for variable-order (e.g., Lagrange) interpolation,
!              the order (1, 2, 3...). Sets the number of 'ghost' zones
!              of data transferred between MPI tasks.
!              scheme
!    iexchtyp: format for exchanging particle data: GPEXCHTYPE_NN=nearest
!    neighbor,
!              suggesting that data is passed between neighboring
!              tasks only; GPEXCHTYPE_VDB form suggests that all tasks retain a
!              copy
!              of the 'master' particle d.b., so that passing between
!              neighbors is not necessary. But the VDB form does require
!              a global reduction, and may be more expensive.
!    iouttup : output type: 0==binary, 1=ASCII
!    bcoll   : if doing binary I/O, do collective (==1); or not (==0)
!    csize   : cache size param for local transposes
!    nstrip  : 'strip-mining' size for local transposes
!    intacc  : compute acceleration internally to class (==1); or not (==0).
!    Storage
!              allocated only if intacc==1.
!    wrtunit : (optional) write particle positions in box units (==1) (i.e.,
!              x,y,z in [0,2.pi]), or in grid units (==0) (x,y,z in [0,N]).
!-----------------------------------------------------------------
    USE var
    USE grid
    USE boxsize
    USE commtypes

    IMPLICIT NONE
    CLASS(GPIC)      ,INTENT(INOUT)     :: this
    INTEGER          ,INTENT   (IN)     :: bcoll,comm,ppc
    INTEGER          ,INTENT   (IN)     :: csize,nstrip,intacc
    INTEGER, OPTIONAL,INTENT   (IN)     :: wrtunit
    INTEGER                             :: disp(3),lens(3),types(3),szreal
    INTEGER          ,INTENT   (IN)     :: iexchtyp,inittype
    INTEGER          ,INTENT   (IN)     :: intorder,iouttyp
    INTEGER                             :: maxparts,i

!    this%icv_      = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)  &
!                     *Dkx*Dky*Dkz/(2*pi)**3
!    this%icv_      = real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)maxparts
    this%icv_      = 1.0_GP/real(ppc,kind=GP)
    maxparts       = ppc*nx*ny*nz

    CALL this%GPart_ctor(comm,maxparts,inittype,0,3,iexchtyp,iouttyp,&
                         bcoll,csize,nstrip,intacc,wrtunit)
! Free memory from unused interpolator
    CALL this%intop_ %ResizeArrays(0,.false.)

!    CALL this%intop_%GPSplineInt_dtor()
    this%intorder_ = intorder
!    CALL this%gpcomm_%GPartComm_ctor(GPCOMM_INTRFC_SF,this%partbuff_,    &
!         this%nd_,this%intorder_/2+1,this%comm_,this%htimers_(GPTIME_COMM))
!    CALL this%gpcomm_%SetCacheParam(csize,nstrip)
!    CALL this%gpcomm_%Init()
!    CALL this%gpcomm_%GFieldComm_ctor()
    CALL this%gpcomm_%AllocRetArrays()

    CALL this%picspl_%GPICSplineInt_ctor(3,this%nd_,this%libnds_,this%lxbnds_, &
         this%tibnds_,this%intorder_,this%intorder_/2+1,this%partbuff_,        &
         this%gpcomm_,this%htimers_(GPTIME_DATAEX),this%htimers_(GPTIME_TRANSP))

    ALLOCATE ( this%prop_  (this%partbuff_) )
    ALLOCATE ( this%weight_(this%partbuff_) )

  END SUBROUTINE GPIC_ctor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_Init(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Init
!  DESCRIPTION: Initializes particle locations before integration.
!               Call after construction.
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    IMPLICIT NONE
    CLASS(GPIC)   ,INTENT(INOUT)    :: this
    INTEGER                         :: j

    IF      ( this%inittype_ .EQ. GPINIT_RANDLOC ) THEN
      CALL GPIC_InitRandSeed (this)
    ELSE IF ( this%inittype_ .EQ. GPINIT_USERLOC ) THEN
      CALL GPIC_InitUserSeed (this)
    ENDIF

  END SUBROUTINE GPIC_Init
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_InitUserSeed(this)
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
    CLASS(GPIC) ,INTENT(INOUT)        :: this

    INTEGER                           :: navg,nl,nowned,nt,j
    INTEGER,ALLOCATABLE,DIMENSION(:)  :: iproc,ilproc
    REAL(KIND=GP)                     :: x,y,z,w
    DOUBLE PRECISION                  :: wm

    ! Note: each record (line) consists of x y z real positions
    ! within [0,NX-1]x[0,NY-1]x[0,NZ-1] box or the equivalent
    ! in box units depending on wrtunit_ class options.
    OPEN(UNIT=5,FILE=trim(this%seedfile_),STATUS='OLD',ACTION='READ',&
         IOSTAT=this%ierr_, IOMSG=this%serr_);
    IF ( this%ierr_ .NE. 0 ) THEN
      WRITE(*,*)'GPIC::InitUserSeed: file:',this%seedfile_,' err: ', trim(this%serr_) 
      STOP
    ENDIF
    READ(5,*,IOSTAT=this%ierr_) nt
    IF ( this%myrank_.eq.0 .AND. nt.NE.this%maxparts_ ) THEN
      WRITE(*,*) 'GPIC_InitUserSeed: Inconsistent seed file: required no. part.=', &
      this%maxparts_,' total listed: ',nt,' file:',this%seedfile_
      STOP
    ENDIF
    READ(5,*,IOSTAT=this%ierr_) x

    nt = 0     ! global part. record counter
    nl = 0     ! local particle counter
    wm = 0.0D0 ! mean particle weight
    DO WHILE ( this%ierr_.EQ.0 )
      READ(5,*,IOSTAT=this%ierr_) x, y, z, w
      IF ( this%ierr_ .NE. 0 ) THEN
!!      WRITE(*,*) 'GPIC::InitUserSeed: terminating read; nt=', nt, ' ierr=',this%ierr_
        EXIT
      ENDIF
      wm = wm + w
      IF ( this%wrtunit_ .EQ. 1 ) THEN ! rescale coordinates to grid units
        x = x*this%invdel_(1)
        y = y*this%invdel_(2)
        z = z*this%invdel_(3)
      ENDIF
      IF ( z.GE.this%lxbnds_(3,1) .AND. z.LT.this%lxbnds_(3,2) .AND. &
           y.GE.this%lxbnds_(2,1) .AND. y.LT.this%lxbnds_(2,2) .AND. &
           x.GE.this%lxbnds_(1,1) .AND. x.LT.this%lxbnds_(1,2) ) THEN
        nl = nl + 1
        this%id_(nl) = nt
        this%px_(nl) = x
        this%py_(nl) = y
        this%pz_(nl) = z
        this%weight_(nl) = w
      ENDIF
      nt = nt + 1
    ENDDO
    CLOSE(5)

    this%nparts_ = nl;
    CALL MPI_ALLREDUCE(nl,nt,1,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
    wm = wm/nt
    wm = wm*this%icv_
    DO j = 1,nl  ! Reescale to unit density
      this%weight_(j) = this%weight_(j)*wm
    END DO
    IF ( this%myrank_.eq.0 .AND. nt.NE.this%maxparts_ ) THEN
      WRITE(*,*) 'GPIC_InitUserSeed: Inconsistent particle count: required no.=', &
      this%maxparts_,' total read: ',nt,' file:',this%seedfile_
      STOP
    ENDIF

    IF (this%iexchtype_.EQ.GPEXCHTYPE_VDB) THEN
      CALL this%gpcomm_%VDBSynch(this%vdb_,this%maxparts_,this%id_, &
                          this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
    END IF

  END SUBROUTINE GPIC_InitUserSeed
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE GPIC_InitLattice(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : InitLattice
!  DESCRIPTION: Initializes particle locations by dividing
!               maxparts evenly among cells, and randomly
!               selecting positions within each cell.
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    USE grid
    USE mpivars

    IMPLICIT NONE
    CLASS(GPIC)   ,INTENT(INOUT)      :: this
    REAL(KIND=GP)                     :: del
    INTEGER                           :: ppc,pps,ib,lag
    INTEGER                           :: i,j,k,ii,jj,kk

    ppc = this%maxparts_/(nx*ny*nz)
    pps = ppc**(1.0/3.0)
    IF (pps*pps*pps .NE. ppc) THEN
      IF ( this%myrank_.eq.0 ) THEN
        WRITE(*,*) 'GPIC_InitLattice: Number of particles per cell &
                                      must be perfect cube'
        STOP
      ENDIF
    END IF
    this%nparts_ = nx*ny*(kend - ksta + 1)*ppc
    ib = nx*ny*(ksta-1)*ppc - 1
    lag = 1
    del = 1.0_GP/pps
    DO i = 1,nx
      DO j = 1,ny
        DO k = ksta,kend
          DO ii = 1,pps
            DO jj = 1,pps
              DO kk = 1,pps
                this%id_(lag) = lag + ib
                this%px_(lag) = (i-1.0_GP) + (ii-0.50_GP)*del
                this%py_(lag) = (j-1.0_GP) + (jj-0.50_GP)*del
                this%pz_(lag) = (k-1.0_GP) + (kk-0.50_GP)*del
                lag = lag + 1
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
    

    CALL this%gpcomm_%VDBSynch(this%vdb_,this%maxparts_,this%id_, &
                          this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
    CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
                          this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
    CALL GPart_GetLocalWrk(this,this%id_,this%px_,this%py_,this%pz_, &
                           this%nparts_,this%vdb_,this%maxparts_)

    IF ( this%wrtunit_ .EQ. 1 ) THEN ! rescale coordinates to box units
       this%ptmp0_(1,:) = this%vdb_(1,:)*this%delta_(1)
       this%ptmp0_(2,:) = this%vdb_(2,:)*this%delta_(2)
       this%ptmp0_(3,:) = this%vdb_(3,:)*this%delta_(3)
       CALL GPart_ascii_write_lag(this,1,'.','xlgInitRndSeed','000',0.0_GP,&
            this%maxparts_,this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:))
    ELSE
       CALL GPart_ascii_write_lag(this,1,'.','xlgInitRndSeed','000',0.0_GP,&
            this%maxparts_,this%vdb_(1,:),this%vdb_(2,:),this%vdb_(3,:))
    ENDIF

    IF ( .NOT.GPart_PartNumConsistent(this,this%nparts_) ) THEN
      IF ( this%myrank_.eq.0 ) THEN
        WRITE(*,*) 'GPIC_InitLattice: Invalid particle after GetLocalWrk call'
        STOP
      ENDIF
    ENDIF

  END SUBROUTINE GPIC_InitLattice
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_InitRandSeed(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : InitRandSeed
!  DESCRIPTION: Initializes particle locations by dividing
!               maxparts evenly among cells, and randomly
!               selecting positions within each cell.
!  ARGUMENTS  :
!    this    : 'this' class instance
!-----------------------------------------------------------------
    USE random
    USE grid
    USE mpivars

    IMPLICIT NONE
    CLASS(GPIC)   ,INTENT(INOUT)      :: this
    REAL(KIND=GP)                     :: r
    INTEGER                           :: ppc,ib,lag,i,j,k,l

    ppc = this%maxparts_/(nx*ny*nz)
    this%nparts_ = nx*ny*(kend - ksta + 1)*ppc
    ib = nx*ny*(ksta-1)*ppc - 1
    lag = 1
    DO i = 1,nx
      DO j = 1,ny
        DO k = ksta,kend
          DO l = 1,ppc
            this%id_(lag) = lag + ib
            CALL prandom_number(r)
            this%px_(lag) = i - 1.0_GP  + r
            CALL prandom_number(r)
            this%py_(lag) = j - 1.0_GP  + r
            CALL prandom_number(r)
            this%pz_(lag) = k - 1.0_GP + r
            lag = lag + 1
          END DO
        END DO
      END DO
    END DO

    CALL this%gpcomm_%VDBSynch(this%vdb_,this%maxparts_,this%id_, &
                          this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
    CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
                          this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
    CALL GPart_GetLocalWrk(this,this%id_,this%px_,this%py_,this%pz_, &
                           this%nparts_,this%vdb_,this%maxparts_)

    IF ( this%wrtunit_ .EQ. 1 ) THEN ! rescale coordinates to box units
       this%ptmp0_(1,:) = this%vdb_(1,:)*this%delta_(1)
       this%ptmp0_(2,:) = this%vdb_(2,:)*this%delta_(2)
       this%ptmp0_(3,:) = this%vdb_(3,:)*this%delta_(3)
       CALL GPart_ascii_write_lag(this,1,'.','xlgInitRndSeed','000',0.0_GP,&
            this%maxparts_,this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:))
    ELSE
       CALL GPart_ascii_write_lag(this,1,'.','xlgInitRndSeed','000',0.0_GP,&
            this%maxparts_,this%vdb_(1,:),this%vdb_(2,:),this%vdb_(3,:))
    ENDIF

    IF ( .NOT.GPart_PartNumConsistent(this,this%nparts_) ) THEN
      IF ( this%myrank_.eq.0 ) THEN
        WRITE(*,*) 'GPIC_InitRandSeed: Invalid particle after GetLocalWrk call'
        STOP
      ENDIF
    ENDIF

  END SUBROUTINE GPIC_InitRandSeed
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_dtor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Explicit destructor for particle-in-cell. Should be called
!  after calling GPart_dtor.
!
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------
    USE var
    USE grid
    USE boxsize
    USE mpivars
    USE commtypes
    USE random

    IMPLICIT NONE
    TYPE(GPIC)       ,INTENT(INOUT)     :: this

    IF ( ALLOCATED ( this%prop_   ) ) DEALLOCATE ( this%prop_   )
    IF ( ALLOCATED ( this%weight_ ) ) DEALLOCATE ( this%weight_ )

  END SUBROUTINE GPIC_dtor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_GetDensity(this,dens)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetDensity
!  DESCRIPTION: Deposits particle density into grid.
!
!  ARGUMENTS  :
!    this     : 'this' class instance
!    dens     : Eulerian field containing particle density (OUT)
!-----------------------------------------------------------------
    USE grid
    USE mpivars

    IMPLICIT NONE
    CLASS(GPIC)  ,INTENT(INOUT)                            :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: dens
    INTEGER                                                :: lag

    DO lag=1,this%nparts_
      this%prop_(lag) = this%weight_(lag)
    END DO

    CALL GPIC_LagToEuler(this,this%prop_,this%nparts_,dens,.true.)

    RETURN

  END SUBROUTINE GPIC_GetDensity
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_EndStageRKK(this,vx,vy,vz,xk)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : EndStageRKK
!  DESCRIPTION: Called at the end of all RK-like stages to
!               complete Lagrangian particle update.

!  ARGUMENTS  :
!    this    : 'this' class instance
!    vz,vy,vz: compoments of velocity field, in real space, partially
!              updated, possibly. These will be overwritten!
!    xk      : multiplicative RK time stage factor
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars
    USE grid

    IMPLICIT NONE
    CLASS(GPIC)  ,INTENT(INOUT)                            :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: vx,vy,vz
    REAL(KIND=GP),INTENT   (IN)                            :: xk
    INTEGER                                                :: j,ng

    ! u(t+dt) = u*: done already

    ! If using nearest-neighbor interface, do particle exchange
    ! between nearest-neighbor tasks BEFORE z-PERIODIZING particle coordinates:
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_NN ) THEN
      IF (this%nprocs_.GT.1) THEN
        CALL GTStart(this%htimers_(GPTIME_COMM))
        CALL this%gpcomm_%PartExchangeV(this%id_,this%px_,this%py_,this%pz_,  &
             this%nparts_,this%lxbnds_(3,1),this%lxbnds_(3,2),GPEXCH_INIT)
        CALL this%gpcomm_%PartExchangeV(this%id_,this%ptmp0_(1,:),            &
             this%ptmp0_(2,:),this%ptmp0_(3,:),this%nparts_,this%lxbnds_(3,1),&
             this%lxbnds_(3,2),GPEXCH_UPDT)
        CALL this%gpcomm_%PartExchangeV(this%id_,this%weight_,this%prop_, &
             this%prop_,this%nparts_,this%lxbnds_(3,1),this%lxbnds_(3,2), &
             GPEXCH_END)
        CALL GTAcc(this%htimers_(GPTIME_COMM))
      END IF
      ! Enforce periodicity in x and y:
      CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,3)
      ! Enforce periodicity in z and ptmp0(3):
      CALL GPart_MakePeriodicZ(this,this%pz_,this%ptmp0_(3,:),this%nparts_)
    ENDIF

    ! If using VDB interface, do synch-up, and get local work:
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN
      ! Enforce periodicity in x, y, & z:
      CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,7)

      IF ( .NOT.GPart_PartNumConsistent(this,this%nparts_) ) THEN
        IF ( this%myrank_.eq.0 ) THEN
          WRITE(*,*) 'GPIC_EndStageRKK: Inconsistent particle count'
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

      CALL this%gpcomm_%VDBSynch(this%ptmp0_,this%maxparts_,this%id_, &
                          this%weight_,this%weight_,this%weight_,this%nparts_,this%ptmp1_)
      CALL GPIC_CopyLocalWrkScalar(this,this%weight_,this%vdb_,this%ptmp0_(1,:),this%maxparts_)

      CALL MPI_ALLREDUCE(this%nparts_,ng,1,MPI_INTEGER,   &
                         MPI_SUM,this%comm_,this%ierr_)

      IF ( this%myrank_.EQ.0 .AND. ng.NE.this%maxparts_) THEN
        WRITE(*,*)'GPIC_EndStageRKK: inconsistent d.b.: expected: ', &
                 this%maxparts_, '; found: ',ng
        CALL GPART_ascii_write_lag(this,1,'.','xlgerr','000',0.0_GP, &
                                   this%maxparts_,this%vdb_)
        STOP
      ENDIF

    ENDIF

    IF ( this%intacc_.EQ.0 ) RETURN

    ! If doing internal acceleration, synch up past time levels:
    CALL GPart_synch_acc(this)

    ! Set t^n+1 velocity based on most recent Lag.particle positions:
    ! NOTE: vx, vy, vz are overwirtten on exit:
    CALL GPIC_EulerToLag(this,this%lvx_,this%nparts_,vx,.true. )
    CALL GPIC_EulerToLag(this,this%lvy_,this%nparts_,vy,.false.)
    CALL GPIC_EulerToLag(this,this%lvz_,this%nparts_,vz,.false.)

!$omp parallel do
    DO j = 1, this%nparts_
      this%vk2_(1,j) = this%lvx_(j)
      this%vk2_(2,j) = this%lvy_(j)
      this%vk2_(3,j) = this%lvz_(j)
    ENDDO

    ! set particle ids owned by task at lag times:
    this%idm_ = 0
    this%npartsm_ = this%nparts_
    this%idm_(1:this%nparts_) = this%id_(1:this%nparts_)

    RETURN

  END SUBROUTINE GPIC_EndStageRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_StepRKK(this, vx, vy, vz, dt, xk)
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
!               Uses picspline method for interpolation.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    vz,vy,vz: compoments of velocity field, in real space, partially
!              updated, possibly. These will be overwritten!
!    dt      : integration timestep
!    xk      : multiplicative RK time stage factor
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars
    USE grid

    IMPLICIT NONE
    CLASS(GPIC)  ,INTENT(INOUT)                            :: this
    INTEGER                                                :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: vx,vy,vz
    REAL(KIND=GP),INTENT   (IN)                            :: dt,xk
    REAL(KIND=GP)                                          :: dtfact
    REAL(KIND=GP),ALLOCATABLE  ,DIMENSION              (:) :: lid,gid

    CALL GTStart(this%htimers_(GPTIME_STEP))

    ! Find F(u*):
    ! ... x:
    CALL GPIC_EulerToLag(this,this%lvx_,this%nparts_,vx,.true.)
    ! ux* <-- ux + dt * F(U*)*xk:
    dtfact = dt*xk*this%invdel_(1)
  !$omp parallel do
    DO j = 1, this%nparts_
      this%px_(j) = this%ptmp0_(1,j) + dtfact*this%lvx_(j)
    ENDDO

    !
    ! ... y:
    ! Exchange bdy data for velocities, so that we
    ! can perform local interpolations:
    CALL GPIC_EulerToLag(this,this%lvy_,this%nparts_,vy,.false.)
    ! uy* <-- uy + dt * F(U*)*xk:
    dtfact = dt*xk*this%invdel_(2)
  !$omp parallel do
    DO j = 1, this%nparts_
      this%py_(j) = this%ptmp0_(2,j) + dtfact*this%lvy_(j)
    ENDDO

    ! ... z:
    ! Exchange bdy data for velocities, so that we
    ! can perform local interpolations:
    CALL GPIC_EulerToLag(this,this%lvz_,this%nparts_,vz,.false.)
    ! uz* <-- uz + dt * F(U*)*xk:
    dtfact = dt*xk*this%invdel_(3)
  !$omp parallel do
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

  ! At this point, the vx, vy, vz should be intact:
    CALL GPIC_EndStageRKK(this,vx,vy,vz,xk)

  END SUBROUTINE GPIC_StepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_EulerToLag(this,lag,nl,evar,doupdate)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : EulerToLag
!  DESCRIPTION: Computes the Eulerian to Lagrangian
!               transformation by interpolating Eulerian field
!               evar to current position of Lagrangian paricles in
!               d.b. Array lag must be large enough to accommodate
!               max. no. particles; no checking is done. Note
!               that 'evar' array must have local dimensions
!               for a real array in GHOST (nx X ny X (kend-ksta+1)).
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
!
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPIC)  ,INTENT(INOUT)                            :: this
    INTEGER      ,INTENT   (IN)                            :: nl
    LOGICAL      ,INTENT   (IN)                            :: doupdate
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: evar
    REAL(KIND=GP),INTENT(INOUT),DIMENSION             (nl) :: lag
    INTEGER                                                :: j

    IF ( doupdate ) THEN
      CALL GTStart(this%htimers_(GPTIME_PUPDATE))
      CALL this%picspl_%PartUpdate(this%px_,this%py_,this%pz_,this%nparts_)
      CALL GTAcc(this%htimers_(GPTIME_PUPDATE))
    ENDIF
    CALL GTStart(this%htimers_(GPTIME_INTERP))
!    CALL GTStart(this%htimers_(GPTIME_SPLINE))
    CALL this%picspl_%CompSpline(evar)
!    CALL GTAcc(this%htimers_(GPTIME_SPLINE))

!    CALL GTStart(this%htimers_(GPTIME_INTERP))
    CALL this%picspl_%DoInterp(lag,nl)
    CALL GTAcc(this%htimers_(GPTIME_INTERP))

  END SUBROUTINE GPIC_EulerToLag
!-----------------------------------------------------------------
!-----------------------------------------------------------------

SUBROUTINE GPIC_LagToEuler(this,lag,nl,evar,doupdate)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : EulerToLag
!  DESCRIPTION: Computes the Lagrangian to Eulerian
!               transformation by depositing Lagrangian particles
!               property in their current position (in d.b.) to
!               Eulerian field evar interpolating Eulerian field
!               Array lag must be large enough to accommodate
!               max. no. particles; no checking is done. Note
!               that 'evar' array must have local dimensions
!               for a real array in GHOST (nx X ny X (kend-ksta+1)).
!
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    lag     : real array with the containing Lagrangian property
!              to be deposited to the grid (IN)
!    nl      : no. Lag. points in lag
!    evar    : Eulerian variable (OUT)
!    doupdate: if true, do interp point update in interpolator; else don't
!
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPIC)  ,INTENT(INOUT)                            :: this
    INTEGER      ,INTENT   (IN)                            :: nl
    LOGICAL      ,INTENT   (IN)                            :: doupdate
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: evar
    REAL(KIND=GP),INTENT(INOUT),DIMENSION             (nl) :: lag
    INTEGER                                                :: j

    IF ( doupdate ) THEN
      CALL GTStart(this%htimers_(GPTIME_PUPDATE))
      CALL this%picspl_%PartUpdate(this%px_,this%py_,this%pz_,this%nparts_)
      CALL GTAcc(this%htimers_(GPTIME_PUPDATE))
    ENDIF

    CALL GTStart(this%htimers_(GPTIME_SPLINE))
    CALL this%picspl_%DoDeposit(lag,nl,evar)
    CALL GTAcc(this%htimers_(GPTIME_SPLINE))

  END SUBROUTINE GPIC_LagToEuler
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_io_write_wgt(this, iunit, dir, spref, nmb, time)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : io_write_wgt
!  DESCRIPTION: Does write of particle velocity d.b. to file. 
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
    CLASS(GPIC)  ,INTENT(INOUT)       :: this
    REAL(KIND=GP),INTENT   (IN)       :: time
    REAL(KIND=GP)                     :: prec(3)
    INTEGER,INTENT(IN)                :: iunit
    INTEGER                           :: fh,j,nt
    INTEGER(kind=MPI_OFFSET_KIND)     :: offset
    CHARACTER(len=*),INTENT(IN)       :: dir
    CHARACTER(len=*),INTENT(IN)       :: nmb
    CHARACTER(len=*),INTENT(IN)       :: spref

!    CALL this%gpcomm_%VDBSynch(this%ptmp0_,this%maxparts_,this%id_, &
!         this%weight_,this%weight_,this%weight_,this%nparts_,this%ptmp1_)
    ! If doing non-collective binary or ascii writes, synch up vector:
    IF ((this%iouttype_.EQ.0 .AND. this%bcollective_.EQ.0).OR.this%iouttype_.EQ.1 ) THEN
      IF (this%iexchtype_.EQ.GPEXCHTYPE_NN) THEN
        CALL this%gpcomm_%LagSynch_t0(this%ptmp0_(1,:),this%maxparts_,this%id_, &
                                      this%weight_,this%nparts_)
      ELSE IF (this%iexchtype_.EQ.GPEXCHTYPE_VDB) THEN
        CALL this%gpcomm_%LagSynch(this%ptmp0_(1,:),this%maxparts_,this%id_,    &
                                   this%weight_,this%nparts_,this%ptmp1_(1,:))
      END IF
    ENDIF

    IF ( this%iouttype_ .EQ. 0 ) THEN
      IF ( this%bcollective_.EQ. 1 ) THEN
        ! pass in the current linear _local_ particle velocity arrays
        CALL GPart_binary_write_lag_co(this,iunit,dir,spref,nmb,time,this%nparts_, &
             this%weight_)
      ELSE
        ! pass in the synched-up VDB (copied to ptmp0_):
        CALL GPart_binary_write_lag_t0(this,iunit,dir,spref,nmb,time,this%maxparts_, &
             this%ptmp0_(1,:))
      ENDIF
    ELSE
      ! pass in the synched-up VDB (copied to ptmp0_):
      CALL GPart_ascii_write_lag(this,iunit,dir,spref,nmb,time,this%maxparts_, &
           this%ptmp0_(1,:))
    ENDIF

  END SUBROUTINE GPIC_io_write_wgt
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_io_read_wgt(this, iunit, dir, spref, nmb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : io_read_wgt
!  DESCRIPTION: Does read of test particle weights from file.
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
    CLASS(GPIC)     ,INTENT(INOUT)            :: this
    REAL(KIND=GP)                             :: time
    INTEGER,INTENT(IN)                        :: iunit
    INTEGER                                   :: ng,j
    CHARACTER(len=*),INTENT   (IN)            :: dir
    CHARACTER(len=*),INTENT   (IN)            :: nmb
    CHARACTER(len=*),INTENT   (IN)            :: spref

    CALL GTStart(this%htimers_(GPTIME_GPREAD))
    IF ( this%iouttype_ .EQ. 0 ) THEN        ! Binary files
      IF ( this%bcollective_ .EQ. 1 ) THEN   ! collective binary
        IF (len_trim(nmb).gt.0 ) THEN
        CALL GPIC_binary_read_pdb_co_scalar(this,iunit, &
        trim(dir) // '/' // trim(spref) // '.' // nmb // '.lag',time,this%ptmp1_(1,:))
        ELSE
!        CALL GPIC_binary_read_pdb_co_scalar(this,iunit, trim(spref),time,this%ptmp1_(1,:))
        CALL GPIC_binary_read_pdb_co_scalar(this,iunit, trim(spref),time,this%ptmp1_(1,:))
        ENDIF
      ELSE                      ! master thread binary
        IF (len_trim(nmb).gt.0 ) THEN
        CALL GPIC_binary_read_pdb_t0_scalar(this,iunit,&
             trim(dir) // '/' // trim(spref) // '.' // nmb // '.lag',time,this%ptmp1_(1,:))
        ELSE
        CALL GPIC_binary_read_pdb_t0_scalar(this,iunit, trim(spref),time,this%ptmp1_(1,:))
        ENDIF
      ENDIF
    ELSE                         ! ASCII files
      IF (len_trim(nmb).gt.0 ) THEN
      CALL GPIC_ascii_read_pdb_scalar(this,iunit,&
           trim(dir) // '/' // trim(spref) // '.' // nmb // '.txt',time,this%ptmp1_(1,:))
      ELSE
      CALL GPIC_ascii_read_pdb_scalar(this,iunit,trim(spref),time,this%ptmp1_(1,:))
      ENDIF
    ENDIF
    CALL GTAcc(this%htimers_(GPTIME_GPREAD))

    IF (this%iexchtype_.EQ.GPEXCHTYPE_VDB) THEN
    ! Store in particle velocity arrays
      CALL GPIC_CopyLocalWrkScalar(this,this%weight_,this%vdb_,this%ptmp1_(1,:),this%maxparts_)
    ELSE IF (this%iexchtype_.EQ.GPEXCHTYPE_NN) THEN
      DO j = 1,this%nparts_
        this%weight_(j) = this%ptmp1_(1,j)
      END DO
    END IF
    CALL MPI_ALLREDUCE(this%nparts_,ng,1,MPI_INTEGER,   &
                       MPI_SUM,this%comm_,this%ierr_)
    IF ( this%myrank_.EQ.0 .AND. ng.NE.this%maxparts_ ) THEN
      WRITE(*,*)'GPIC_io_read_wgt: inconsistent d.b.: expected: ', &
                 this%maxparts_, '; found: ',ng
      STOP
    ENDIF

  END SUBROUTINE GPIC_io_read_wgt
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_binary_read_pdb_co_scalar(this,iunit,sfile,time,pdb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : binary_read_pdb_co_scalar
!  DESCRIPTION: Does read of binary Lagrangian particle scalar data 
!               from file, collectively.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iunit   : unit number
!    sfile   : fully resolved file name
!    pdb     : part. weights
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPIC)  ,INTENT(INOUT)               :: this
    REAL(KIND=GP)                             :: rvar,time
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:)  :: pdb
    INTEGER,INTENT(IN)                        :: iunit
    INTEGER                                   :: fh,i,j,lc,nerr,szreal,nb,nr
    INTEGER(kind=MPI_OFFSET_KIND)             :: offset
    CHARACTER(len=*),INTENT   (IN)            :: sfile

    CALL MPI_TYPE_SIZE(GC_REAL,szreal,this%ierr_)
    CALL MPI_FILE_OPEN(this%comm_,trim(sfile),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,this%ierr_)
    IF ( this%ierr_ .NE. MPI_SUCCESS ) THEN
      CALL MPI_ERROR_STRING(this%ierr_, this%serr_, nerr,ierr);
      WRITE(*,*) 'GPIC_binary_read_pdb_co: Error reading opening : ', trim(sfile),&
                trim(this%serr_)
      STOP
    ENDIF

    ! Must read part. data from correct spot in file:
    offset = 0
    CALL MPI_FILE_READ_AT_ALL(fh,offset,rvar,1,GC_REAL,this%istatus_,this%ierr_)
!  no.parts
    IF ( int(rvar).NE.this%maxparts_ ) THEN
      WRITE(*,*) 'GPIC_binary_read_pdb_co_scalar: Attempt to read incorrect number of particles: required:',&
                  this%maxparts_,' no. read: ',int(rvar)
      WRITE(*,*) 'GPIC_binary_read_pdb_co_scalar: Error reading: ', trim(sfile)
      STOP
    ENDIF
    offset = szreal
    CALL MPI_FILE_READ_AT_ALL(fh,offset,rvar,1,GC_REAL,this%istatus_,this%ierr_) ! time
    offset = 2*szreal
    IF (this%iexchtype_.EQ.GPEXCHTYPE_NN) THEN
      nb = 0
      nr = this%maxparts_/this%nprocs_
      i  = 1
      DO WHILE ((this%ierr_.EQ.MPI_SUCCESS) .AND. (nb.LT.this%maxparts_))
        nr = MIN(nr, this%maxparts_-nb)
        CALL MPI_FILE_READ_AT_ALL(fh,offset,this%ptmp1_(1,:),nr,GC_REAL,this%istatus_,this%ierr_) ! PDB
        offset = offset + nr*szreal
        DO j = 1,nr
          IF ((i.LE.this%nparts_).AND.(this%id_(i).EQ.(j+nb-1))) THEN
            pdb(i) = this%ptmp1_(1,j)
            i = i + 1
          END IF
        END DO
        nb = nb + nr
      END DO
    ELSE IF (this%iexchtype_.EQ.GPEXCHTYPE_VDB) THEN
      CALL MPI_FILE_READ_AT_ALL(fh,offset,pdb,this%maxparts_,GC_REAL,this%istatus_,this%ierr_)! PDB
    END IF
    CALL MPI_FILE_CLOSE(fh,this%ierr_)

  END SUBROUTINE GPIC_binary_read_pdb_co_scalar
!-----------------------------------------------------------------
!-----------------------------------------------------------------

 SUBROUTINE GPIC_binary_read_pdb_t0_scalar(this,iunit,sfile,time,pdb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : binary_read_pdb_t0
!  DESCRIPTION: Does binary read of Lagrangian position d.b. from file
!               only from MPI task 0, and broadcast to all other tasks.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iunit   : unit number
!    sfile   : fully resolved file name
!    time    : real time
!    pdb     : part. d.b. in (3,maxparts) array
!
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPIC)  ,INTENT(INOUT)              :: this
    REAL(KIND=GP),INTENT(INOUT)              :: time
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(:) :: pdb
    REAL(KIND=GP)                            :: fnt
    INTEGER      ,INTENT   (IN)              :: iunit
    INTEGER                                  :: j
    CHARACTER(len=*),INTENT(IN)              :: sfile

    ! Read global VDB, with time header, indexed only
    ! by time index: dir/spref.TTT.lag:
    IF ( this%myrank_.EQ.0 ) THEN
      OPEN(iunit,file=trim(sfile),status='old',access='stream', &
           form='unformatted',iostat=this%ierr_)
      IF ( this%ierr_.NE.0 ) THEN
        WRITE(*,*)'GPIC_binary_read_pdb_t0_scalar: could not open file for reading:',&
        trim(sfile)
        STOP
      ENDIF

      REWIND(iunit)
      READ(iunit) fnt
      READ(iunit) time
      IF ( int(fnt).NE.this%maxparts_ ) THEN
        WRITE(*,*)this%myrank_, &
          ': GPIC_binary_read_pdb_t0_scalar: particle inconsistency: no. required=',&
          this%maxparts_,' no. found=',int(fnt), &
          ' file=',trim(sfile)
        STOP
      ENDIF
      READ(iunit) pdb
      CLOSE(iunit)
    ENDIF
 
    IF (this%iexchtype_.EQ.GPEXCHTYPE_VDB) THEN
      CALL MPI_BCAST(pdb,this%maxparts_,GC_REAL,0,this%comm_,this%ierr_)
      IF ( this%ierr_.NE.MPI_SUCCESS ) THEN
        WRITE(*,*)this%myrank_, ': GPIC_binary_read_pdb_t0_scalar: Broadcast failed: file=',&
        trim(sfile)
      ENDIF
    ELSE IF (this%iexchtype_.EQ.GPEXCHTYPE_NN) THEN
      CALL this%gpcomm_%PartScatterV(this%id_,pdb,pdb,pdb,this%nparts_,this%tmpint_)
    END IF

  END SUBROUTINE GPIC_binary_read_pdb_t0_scalar
!-----------------------------------------------------------------
!-----------------------------------------------------------------

 SUBROUTINE GPIC_ascii_read_pdb_scalar(this,iunit,sfile,time,pdb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : ascii_read_pdb
!  DESCRIPTION: Does ASCII read of Lagrangian position d.b. from file.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    iunit   : unit number
!    sfile   : fully resolved file name
!    time    : real time
!    pdb     : part. d.b. in (4,maxparts) array
!
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(GPIC)  ,INTENT(INOUT)       :: this
    REAL(KIND=GP),INTENT(INOUT)       :: time
    REAL(KIND=GP),INTENT(INOUT)       :: pdb(this%maxparts_)
    INTEGER      ,INTENT   (IN)       :: iunit
    INTEGER                           :: j,nt
    CHARACTER(len=*),INTENT(IN)       :: sfile

    ! Read global VDB, with time header, indexed only
    ! by time index: dir/spref.TTT.txt:
    IF ( this%myrank_.EQ.0 ) THEN
      OPEN(iunit,file=trim(sfile),status='old',form='formatted',iostat=this%ierr_)
      IF ( this%ierr_.NE.0 ) THEN
        WRITE(*,*)'GPIC_ascii_read_pdb_scalar: could not open file for reading: ',&
        trim(sfile)
        STOP
      ENDIF
      READ(iunit,*,iostat=this%ierr_) nt
      READ(iunit,*,iostat=this%ierr_) time
      IF ( nt.LT.this%maxparts_ ) THEN
        WRITE(*,*)this%myrank_, &
          ': GPIC_ascii_read_pdb_scalar: particle inconsistency: no. required=',&
          this%maxparts_,' no. found=',nt, &
          ' file=',trim(sfile)
        STOP
      ENDIF
      DO j = 1, this%maxparts_
        READ(iunit,*,iostat=this%ierr_) pdb(j)
  600   FORMAT(3(E23.15,1X))
      ENDDO
      CLOSE(iunit)
    ENDIF
    CALL MPI_BCAST(pdb,this%maxparts_,GC_REAL,0,this%comm_,this%ierr_)
    IF ( this%ierr_.NE.MPI_SUCCESS ) THEN
        WRITE(*,*)this%myrank_, ': GPIC_ascii_read_pdb_scalar: Broadcast failed:file=',&
        trim(sfile)
    ENDIF

  END SUBROUTINE GPIC_ascii_read_pdb_scalar
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_CopyLocalWrkScalar(this,l,gvdb,vgvdb,ngvdb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GPIC_CopyLocalWrkScalar
!  DESCRIPTION: Updates records of the VDB.
!  ARGUMENTS  :
!    this    : 'this' class instance (IN)
!    l       : local part. d.b. vector
!    gvdb    : global VDB containing part. position records. Location
!              gives particle id.
!    vgvdb   : global VDB containing part. property records (scalar)
!    ngvdb   : no. records in global VDB
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes

    IMPLICIT NONE
    CLASS(GPIC)  ,INTENT(INOUT)                           :: this
    INTEGER      ,INTENT   (IN)                           :: ngvdb
    INTEGER                                               :: i,j,nll
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(this%maxparts_) :: l
    REAL(KIND=GP),INTENT   (IN),DIMENSION(3,ngvdb)        :: gvdb
    REAL(KIND=GP),INTENT   (IN),DIMENSION(ngvdb)          :: vgvdb

    nll = 0
    DO j = 1, ngvdb
      IF ( gvdb(3,j).GE.this%lxbnds_(3,1) .AND. gvdb(3,j).LT.this%lxbnds_(3,2) ) THEN 
        nll = nll + 1
        l (nll) = vgvdb(j)
      ENDIF
    ENDDO

  END SUBROUTINE GPIC_CopyLocalWrkScalar
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE GPIC_ResizeArrays(this,new_size,onlyinc,exc)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Resize_Arrays
!  DESCRIPTION: Resize all arrays in the VGPIC class (including 
!               subclases, i.e. communicator, spline)
!  ARGUMENTS  :
!    this    : 'this' class instance
!    new_size: new number of particles
!    onlyinc : if true, will only resize to increase array size
!-----------------------------------------------------------------
!$  USE threads
 
    IMPLICIT NONE
    CLASS(GPIC) ,INTENT(INOUT)                         :: this
    INTEGER     ,INTENT(IN)                            :: new_size
    LOGICAL     ,INTENT(IN)                            :: onlyinc
    LOGICAL     ,INTENT(IN)   ,OPTIONAL                :: exc
    INTEGER                                            :: n

    CALL GPart_ResizeArrays(this,new_size,onlyinc)

    n = SIZE(this%prop_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank1(this%prop_,new_size,.false.)
    END IF
    n = SIZE(this%weight_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank1(this%weight_,new_size,.true.)
    END IF
 
    IF (PRESENT(exc)) THEN
      IF (exc) RETURN    ! Skip subclass resizing
    END IF

! Field communicator resized within picspline
    CALL this%picspl_%ResizeArrays(new_size,onlyinc)

    RETURN

  END SUBROUTINE GPIC_ResizeArrays
!-----------------------------------------------------------------
!-----------------------------------------------------------------


!=================================================================
! VGPIC SUBROUTINES
!=================================================================

   SUBROUTINE VGPIC_InitUserSeed(this)
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
    CLASS(VGPIC) ,INTENT(INOUT)       :: this

    INTEGER                           :: navg,nl,nowned,nt,j
    INTEGER,ALLOCATABLE,DIMENSION(:)  :: iproc,ilproc
    REAL(KIND=GP)                     :: x,y,z,w,vx,vy,vz
    DOUBLE PRECISION                  :: wm

    ! Note: each record (line) consists of x y z w vx vy vz
    ! The x y z are real positions within 
    ! [0,NX-1]x[0,NY-1]x[0,NZ-1] box or the equivalent
    ! in box units depending on wrtunit_ class options.
    ! The w is a real particle weight, which will be
    ! reescaled to match unit mean particle density.
    ! The vx vy vz are real velocities, using the same
    ! units as the x y z positions
    OPEN(UNIT=5,FILE=trim(this%seedfile_),STATUS='OLD',ACTION='READ',&
         IOSTAT=this%ierr_, IOMSG=this%serr_);
    IF ( this%ierr_ .NE. 0 ) THEN
      WRITE(*,*)'VGPIC::InitUserSeed: file:',this%seedfile_,' err: ', trim(this%serr_) 
      STOP
    ENDIF
    READ(5,*,IOSTAT=this%ierr_) nt
    IF ( this%myrank_.eq.0 .AND. nt.NE.this%maxparts_ ) THEN
      WRITE(*,*) 'VGPIC_InitUserSeed: Inconsistent seed file: required no. part.=', &
      this%maxparts_,' total listed: ',nt,' file:',this%seedfile_
      STOP
    ENDIF
    READ(5,*,IOSTAT=this%ierr_) x

    nt = 0     ! global part. record counter
    nl = 0     ! local particle counter
    wm = 0.0D0 ! mean particle weight
    DO WHILE ( this%ierr_.EQ.0 )
      READ(5,*,IOSTAT=this%ierr_) x, y, z, w, vx, vy, vz
      IF ( this%ierr_ .NE. 0 ) THEN
      WRITE(*,*) 'VGPIC::InitUserSeed: terminating read; nt=', nt, ' ierr=',this%ierr_
        EXIT
      ENDIF
      wm = wm + w
      IF ( this%wrtunit_ .EQ. 1 ) THEN ! rescale coordinates to grid units
        x  =  x*this%invdel_(1)
        y  =  y*this%invdel_(2)
        z  =  z*this%invdel_(3)
      ELSE
        vx = vx*this%delta_(1)
        vy = vy*this%delta_(2)
        vz = vz*this%delta_(3)
      ENDIF
      IF ( z.GE.this%lxbnds_(3,1) .AND. z.LT.this%lxbnds_(3,2) .AND. &
           y.GE.this%lxbnds_(2,1) .AND. y.LT.this%lxbnds_(2,2) .AND. &
           x.GE.this%lxbnds_(1,1) .AND. x.LT.this%lxbnds_(1,2) ) THEN
        nl = nl + 1
        this%id_(nl)  = nt
        this%px_(nl)  = x
        this%py_(nl)  = y
        this%pz_(nl)  = z
        this%pvx_(nl) = vx
        this%pvy_(nl) = vy
        this%pvz_(nl) = vz
        this%weight_(nl) = w
      ENDIF
      nt = nt + 1
    ENDDO
    CLOSE(5)

    this%nparts_ = nl;
    CALL MPI_ALLREDUCE(nl,nt,1,MPI_INTEGER,MPI_SUM,this%comm_,this%ierr_)
    wm = wm/nt
    wm = wm*this%icv_
    DO j = 1,nl  ! Reescale to unit density
      this%weight_(j) = this%weight_(j)*wm
    END DO
    IF ( this%myrank_.eq.0 .AND. nt.NE.this%maxparts_ ) THEN
      WRITE(*,*) 'VGPIC_InitUserSeed: Inconsistent particle count: required no.=', &
      this%maxparts_,' total read: ',nt,' file:',this%seedfile_
      STOP
    ENDIF

    IF (this%iexchtype_.EQ.GPEXCHTYPE_VDB) THEN
      CALL this%gpcomm_%VDBSynch(this%vdb_,this%maxparts_,this%id_, &
                          this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
    END IF

  END SUBROUTINE VGPIC_InitUserSeed
!-----------------------------------------------------------------
!-----------------------------------------------------------------


  SUBROUTINE VGPIC_io_write_pdbv(this, iunit, dir, spref, nmb, time)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : io_write_pdbv
!  DESCRIPTION: Does write of particle velocity d.b. to file.
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
    CLASS(VGPIC),INTENT(INOUT)        :: this
    REAL(KIND=GP),INTENT   (IN)       :: time
    REAL(KIND=GP)                     :: prec(3)
    INTEGER,INTENT(IN)                :: iunit
    INTEGER                           :: fh,j,nt
    INTEGER(kind=MPI_OFFSET_KIND)     :: offset
    CHARACTER(len=*),INTENT(IN)       :: dir
    CHARACTER(len=*),INTENT(IN)       :: nmb
    CHARACTER(len=*),INTENT(IN)       :: spref

!    CALL this%gpcomm_%VDBSynch(this%ptmp0_,this%maxparts_,this%id_, &
!         this%pvx_,this%pvy_,this%pvz_,this%nparts_,this%ptmp1_)

    ! If doing non-collective binary or ascii writes, synch up vector:
    IF ((this%iouttype_.EQ.0 .AND.this%bcollective_.EQ.0).OR.this%iouttype_.EQ.1 ) THEN
      IF (this%iexchtype_.EQ.GPEXCHTYPE_NN) THEN
        CALL this%gpcomm_%VDBSynch_t0(this%ptmp0_,this%maxparts_,this%id_, &
                                      this%pvx_,this%pvy_,this%pvz_,this%nparts_)
      ELSE IF (this%iexchtype_.EQ.GPEXCHTYPE_VDB) THEN
        CALL this%gpcomm_%VDBSynch(this%ptmp0_,this%maxparts_,this%id_,     &
                                   this%pvx_,this%pvy_,this%pvz_,this%nparts_,this%ptmp1_)
      END IF
    ENDIF

    IF ( this%iouttype_ .EQ. 0 ) THEN
      IF ( this%bcollective_.EQ. 1 ) THEN
        ! pass in the current linear _local_ particle velocity arrays
        CALL GPart_binary_write_lag_co(this,iunit,dir,spref,nmb,time,this%nparts_, &
             this%pvx_,this%pvy_,this%pvz_)
      ELSE
        ! pass in the synched-up VDB (copied to ptmp0_):
        CALL GPart_binary_write_lag_t0(this,iunit,dir,spref,nmb,time,this%maxparts_, &
             this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:))
      ENDIF
    ELSE
      ! pass in the synched-up VDB (copied to ptmp0_):
      CALL GPart_ascii_write_lag(this,iunit,dir,spref,nmb,time,this%maxparts_, &
           this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:))
    ENDIF

  END SUBROUTINE VGPIC_io_write_pdbv
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE VGPIC_io_read_pdbv(this, iunit, dir, spref, nmb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : io_read_pdbv
!  DESCRIPTION: Does read of test particle velocity from file.
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
    CLASS(VGPIC)    ,INTENT(INOUT)            :: this
    REAL(KIND=GP)                             :: time
    INTEGER,INTENT(IN)                        :: iunit
    INTEGER                                   :: ng,j
    CHARACTER(len=*),INTENT   (IN)            :: dir
    CHARACTER(len=*),INTENT   (IN)            :: nmb
    CHARACTER(len=*),INTENT   (IN)            :: spref

    CALL GTStart(this%htimers_(GPTIME_GPREAD))
    IF ( this%iouttype_ .EQ. 0 ) THEN        ! Binary files
      IF ( this%bcollective_ .EQ. 1 ) THEN   ! collective binary
        IF (len_trim(nmb).gt.0 ) THEN
          CALL GPart_binary_read_pdb_co(this,iunit, &
               trim(dir) // '/' // trim(spref) // '.' // nmb // '.lag', &
               time,this%ptmp0_)
        ELSE
        CALL GPart_binary_read_pdb_co(this,iunit,trim(spref),time,this%ptmp0_)
        ENDIF
      ELSE                      ! master thread binary
        IF (len_trim(nmb).gt.0 ) THEN
          CALL GPart_binary_read_pdb_t0(this,iunit,&
            trim(dir) // '/' // trim(spref) // '.' // nmb // '.lag',time,this%ptmp0_)
        ELSE
        CALL GPart_binary_read_pdb_t0(this,iunit, trim(spref),time,this%ptmp0_)
        ENDIF
      ENDIF
    ELSE                         ! ASCII files
      IF (len_trim(nmb).gt.0 ) THEN
        CALL GPart_ascii_read_pdb (this,iunit,&
            trim(dir) // '/' // trim(spref) // '.' // nmb // '.txt',time,this%ptmp0_)
      ELSE
      CALL GPart_ascii_read_pdb (this,iunit,trim(spref),time,this%ptmp0_)
      ENDIF
    ENDIF
    CALL GTAcc(this%htimers_(GPTIME_GPREAD))

    IF (this%iexchtype_.EQ.GPEXCHTYPE_VDB) THEN
    ! Store in particle velocity arrays
      CALL GPart_CopyLocalWrk(this,this%pvx_,this%pvy_,this%pvz_, &
                              this%vdb_,this%ptmp0_,this%maxparts_)
    ELSE
      DO j = 1,this%nparts_
        this%pvx_(j) = this%ptmp0_(1,j)
        this%pvy_(j) = this%ptmp0_(2,j)
        this%pvz_(j) = this%ptmp0_(3,j)
      END DO
    END IF

    CALL MPI_ALLREDUCE(this%nparts_,ng,1,MPI_INTEGER,   &
                       MPI_SUM,this%comm_,this%ierr_)
    IF ( this%myrank_.EQ.0 .AND. ng.NE.this%maxparts_ ) THEN
      WRITE(*,*)'VGPIC_io_read_pdbv: inconsistent d.b.: expected: ', &
                 this%maxparts_, '; found: ',ng
      STOP
    ENDIF

  END SUBROUTINE VGPIC_io_read_pdbv
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE VGPIC_GetFlux(this,jx,jy,jz)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetFlux
!  DESCRIPTION: Deposits particle current density into grid.
!
!  ARGUMENTS  :
!    this     : 'this' class instance
!    jx,jy,jz : Eulerian field containing particle current
!               density (OUT)
!-----------------------------------------------------------------
    USE grid
    USE mpivars

    IMPLICIT NONE
    CLASS(VGPIC)  ,INTENT(INOUT)                           :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: jx,jy,jz
    INTEGER                                                :: lag

! x-coord
    DO lag=1,this%nparts_
      this%prop_(lag) = this%weight_(lag)*this%pvx_(lag)
    END DO

    CALL GPIC_LagToEuler(this,this%prop_,this%nparts_,jx,.true.)

! y-coord
    DO lag=1,this%nparts_
      this%prop_(lag) = this%weight_(lag)*this%pvy_(lag)
    END DO
    CALL GPIC_LagToEuler(this,this%prop_,this%nparts_,jy,.false.)

! z-coord
    DO lag=1,this%nparts_
      this%prop_(lag) = this%weight_(lag)*this%pvz_(lag)
    END DO
    CALL GPIC_LagToEuler(this,this%prop_,this%nparts_,jz,.false.)

    RETURN

  END SUBROUTINE VGPIC_GetFlux
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE VGPIC_GetMoment(this,field,ox,oy,oz)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetMoment
!  DESCRIPTION: Deposits particle moment density of chosen component
!
!  ARGUMENTS  :
!    this     : 'this' class instance
!    field    : Eulerian field containing particle moment
!               density (OUT)
!    ox,oy,oz : Order of the moment to calculate in each component
!-----------------------------------------------------------------
    USE grid
    USE mpivars

    IMPLICIT NONE
    CLASS(VGPIC) ,INTENT(INOUT)                            :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: field
    INTEGER      ,INTENT(IN)                               :: ox,oy,oz
    INTEGER                                                :: lag

    DO lag=1,this%nparts_
      this%prop_(lag) = this%weight_(lag)*this%pvx_(lag)**ox &
                                         *this%pvy_(lag)**oy &
                                         *this%pvz_(lag)**oz
    END DO

    CALL GPIC_LagToEuler(this,this%prop_,this%nparts_,field,.false.)

    RETURN

  END SUBROUTINE VGPIC_GetMoment
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE VGPIC_ResizeArrays(this,new_size,onlyinc,exc)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Resize_Arrays
!  DESCRIPTION: Resize all arrays in the VGPIC class (including 
!               subclases, i.e. communicator, spline)
!  ARGUMENTS  :
!    this    : 'this' class instance
!    new_size: new number of particles
!    onlyinc : if true, will only resize to increase array size
!-----------------------------------------------------------------
!$  USE threads
 
    IMPLICIT NONE
    CLASS(VGPIC) ,INTENT(INOUT)                         :: this
    INTEGER      ,INTENT(IN)                            :: new_size
    LOGICAL      ,INTENT(IN)                            :: onlyinc
    LOGICAL      ,INTENT(IN)   ,OPTIONAL                :: exc
    INTEGER                                             :: n

    CALL GPIC_ResizeArrays(this,new_size,onlyinc,exc)

    n = SIZE(this%pvx_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank1(this%pvx_,new_size,.true.)
    END IF
    n = SIZE(this%pvy_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank1(this%pvy_,new_size,.true.)
    END IF
    n = SIZE(this%pvz_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank1(this%pvz_,new_size,.true.)
    END IF

    n = SIZE(this%ttmp0_,2)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank2(this%ttmp0_,new_size,.true.)
    END IF

    RETURN

  END SUBROUTINE VGPIC_ResizeArrays
!-----------------------------------------------------------------
!-----------------------------------------------------------------


!=================================================================
! ChargPIC SUBROUTINES
!=================================================================

  SUBROUTINE ChargPIC_ctor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Explicit constructor for test particles. Should be called
!  after calling GPart_ctor. Calling GPIC_ctor is not required.
!
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------
    USE var
    USE grid
    USE boxsize
    USE commtypes

    IMPLICIT NONE
    CLASS(ChargPIC)   , INTENT(INOUT)   :: this

    ALLOCATE(this%pvx_     (this%partbuff_))
    ALLOCATE(this%pvy_     (this%partbuff_))
    ALLOCATE(this%pvz_     (this%partbuff_))
    ALLOCATE(this%lbx_     (this%partbuff_))
    ALLOCATE(this%lby_     (this%partbuff_))
    ALLOCATE(this%lbz_     (this%partbuff_))
    ALLOCATE(this%lfx_     (this%partbuff_))
    ALLOCATE(this%lfy_     (this%partbuff_))
    ALLOCATE(this%lfz_     (this%partbuff_))
    ALLOCATE(this%ttmp0_ (3,this%partbuff_))

  END SUBROUTINE ChargPIC_ctor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE ChargPIC_dtor(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Explicit destructor for test particles. Should be called
!  after GPart_dtor.
!
!  ARGUMENTS:
!    this    : 'this' class instance
!-----------------------------------------------------------------

    IMPLICIT NONE
    TYPE(ChargPIC)    ,INTENT(INOUT)    :: this

    IF ( ALLOCATED  (this%prop_) ) DEALLOCATE  (this%prop_)
    IF ( ALLOCATED   (this%pvx_) ) DEALLOCATE   (this%pvx_)
    IF ( ALLOCATED   (this%pvy_) ) DEALLOCATE   (this%pvy_)
    IF ( ALLOCATED   (this%pvz_) ) DEALLOCATE   (this%pvz_)
    IF ( ALLOCATED   (this%lbx_) ) DEALLOCATE   (this%lbx_)
    IF ( ALLOCATED   (this%lby_) ) DEALLOCATE   (this%lby_)
    IF ( ALLOCATED   (this%lbz_) ) DEALLOCATE   (this%lbz_)
    IF ( ALLOCATED   (this%lfx_) ) DEALLOCATE   (this%lfx_)
    IF ( ALLOCATED   (this%lfy_) ) DEALLOCATE   (this%lfy_)
    IF ( ALLOCATED   (this%lfz_) ) DEALLOCATE   (this%lfz_)
    IF ( ALLOCATED (this%ttmp0_) ) DEALLOCATE (this%ttmp0_)

  END SUBROUTINE ChargPIC_dtor
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE ChargPIC_SetStepRKK(this)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : SetStepRKK
!  DESCRIPTION: Initializes an explicit integration timestep for
!               the particle velocity. Must be called at the start
!               of the RK stage execution, together with
!               GPart_SetStepRKK.
!
!  ARGUMENTS  :
!    this     : 'this' class instance
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(ChargPIC)   , INTENT(INOUT) :: this

    INTEGER                           :: j

    ! Initialize solution, v (particle velocity):
    ! v* <-- v:

    ! Cycle over JST loop to update state:
!$omp parallel do
    DO j = 1, this%nparts_
       this%ttmp0_(1,j) = this%pvx_(j)  ! ux_0
       this%ttmp0_(2,j) = this%pvy_(j)  ! uy_0
       this%ttmp0_(3,j) = this%pvz_(j)  ! uz_0
    ENDDO

  END SUBROUTINE ChargPIC_SetStepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE ChargPIC_StepBoris(this, Ex, Ey, Ez, Bx, By, Bz, dt, o)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Step_chargedpicBoris
!  DESCRIPTION: Carries out one stage of explicit RK-like time
!               integration step.  Intended for explicit step within
!               an outer stepper method of the form:
!
!               X = X_0 + dt * V[X(t),t] * xk,
!               V = V_0 + dt * F[V(X(t)),E(X(t)),B(X(t))] * xk,
!
!               where F is the electromagnetic force on the particle.
!               Note that the vx, vy, vz, will be overwritten here.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    Ez,Ey,Ez: compoments of electric field, in real space, partially
!              updated, possibly. These will be overwritten!
!    Bz,By,Bz: compoments of magnetic field in real space
!    dt      : integration timestep
!    xk      : multiplicative RK time stage factor
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(ChargPIC)    ,INTENT(INOUT)                      :: this
    INTEGER                                                :: i,j
    INTEGER      ,INTENT   (IN)                            :: o
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: Bx,By,Bz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: Ex,Ey,Ez
    REAL(KIND=GP),INTENT   (IN)                            :: dt
    REAL(KIND=GP)                                          :: dtv,dtx,dty,dtz
    REAL(KIND=GP)                                          :: fact,inprod,b2
    REAL(KIND=GP)                                          :: c1,c2,c3
    REAL(KIND=GP), ALLOCATABLE, DIMENSION              (:) :: lid,gid

    dtv = dt*0.50_GP/real(o,kind=GP)
    dtx = dtv*this%invdel_(1)
    dty = dtv*this%invdel_(2)
    dtz = dtv*this%invdel_(3)
    CALL GTStart(this%htimers_(GPTIME_STEP))

    ! Find E and B:
    CALL GPIC_EulerToLag(this,this%lfx_,this%nparts_,Ex,.true. )
    CALL GPIC_EulerToLag(this,this%lfy_,this%nparts_,Ey,.false.)
    CALL GPIC_EulerToLag(this,this%lfz_,this%nparts_,Ez,.false.)
    CALL GPIC_EulerToLag(this,this%lbx_,this%nparts_,Bx,.false.)
    CALL GPIC_EulerToLag(this,this%lby_,this%nparts_,By,.false.)
    CALL GPIC_EulerToLag(this,this%lbz_,this%nparts_,Bz,.false.)
!$omp parallel do if(this%nparts_.ge.NMIN_OMP)
    DO j = 1, this%nparts_
       this%lbx_(j) = dtv*this%lbx_(j)
       this%lby_(j) = dtv*this%lby_(j)
       this%lbz_(j) = dtv*this%lbz_(j)

       this%pvx_(j) = this%ttmp0_(1,j) + dtv*this%lfx_(j)
       this%pvy_(j) = this%ttmp0_(2,j) + dtv*this%lfy_(j)
       this%pvz_(j) = this%ttmp0_(3,j) + dtv*this%lfz_(j)

       b2     = this%lbx_(j)*this%lbx_(j)+this%lby_(j)*this%lby_(j)+this%lbz_(j)*this%lbz_(j)
       fact   = 2.0_GP/(1.0_GP + b2)
       inprod = this%pvx_(j)*this%lbx_(j)+this%pvy_(j)*this%lby_(j)+this%pvz_(j)*this%lbz_(j)
       c1 = this%pvy_(j)*this%lbz_(j) - this%pvz_(j)*this%lby_(j)  !(vxb)_x
       c2 = this%pvz_(j)*this%lbx_(j) - this%pvx_(j)*this%lbz_(j)  !(vxb)_y
       c3 = this%pvx_(j)*this%lby_(j) - this%pvy_(j)*this%lbx_(j)  !(vxb)_z
       this%pvx_(j) = this%pvx_(j) + fact*(inprod*this%lbx_(j) - this%pvx_(j)*b2 + c1)
       this%pvy_(j) = this%pvy_(j) + fact*(inprod*this%lby_(j) - this%pvy_(j)*b2 + c2)
       this%pvz_(j) = this%pvz_(j) + fact*(inprod*this%lbz_(j) - this%pvz_(j)*b2 + c3)

       this%pvx_(j) = this%pvx_(j) + dtv*this%lfx_(j)
       this%pvy_(j) = this%pvy_(j) + dtv*this%lfy_(j)
       this%pvz_(j) = this%pvz_(j) + dtv*this%lfz_(j)

       this%px_(j) = this%ptmp0_(1,j) + dtx*(this%pvx_(j) + this%ttmp0_(1,j))
       this%py_(j) = this%ptmp0_(2,j) + dty*(this%pvy_(j) + this%ttmp0_(2,j))
       this%pz_(j) = this%ptmp0_(3,j) + dtz*(this%pvz_(j) + this%ttmp0_(3,j))
    ENDDO

    CALL GTAcc(this%htimers_(GPTIME_STEP))

    CALL ChargPIC_EndStageRKK(this,Ex,Ey,Ez,1.0_GP)

  END SUBROUTINE ChargPIC_StepBoris
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE ChargPIC_StepRKK(this, Ex, Ey, Ez, Bx, By, Bz, dt, xk)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Step_chargedpic
!  DESCRIPTION: Carries out one stage of explicit RK-like time
!               integration step.  Intended for explicit step within
!               an outer stepper method of the form:
!
!               X = X_0 + dt * V[X(t),t] * xk,
!               V = V_0 + dt * F[V(X(t)),E(X(t)),B(X(t))] * xk,
!
!               where F is the electromagnetic force on the particle.
!               Note that the vx, vy, vz, will be overwritten here.
!  ARGUMENTS  :
!    this    : 'this' class instance
!    Ez,Ey,Ez: compoments of electric field, in real space, partially
!              updated, possibly. These will be overwritten!
!    Bz,By,Bz: compoments of magnetic field in real space
!    dt      : integration timestep
!    xk      : multiplicative RK time stage factor
!-----------------------------------------------------------------
    USE grid
    USE fprecision
    USE commtypes
    USE mpivars

    IMPLICIT NONE
    CLASS(ChargPIC)    ,INTENT(INOUT)                      :: this
    INTEGER                                                :: i,j
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: Bx,By,Bz
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: Ex,Ey,Ez
    REAL(KIND=GP),INTENT   (IN)                            :: dt,xk
    REAL(KIND=GP)                                          :: dtfact
    REAL(KIND=GP)                                          :: dtv
    REAL(KIND=GP), ALLOCATABLE, DIMENSION              (:) :: lid,gid

    dtv    = dt*xk
    CALL GTStart(this%htimers_(GPTIME_STEP))

    ! Find F(u*):
    CALL GPIC_EulerToLag(this,this%lbx_,this%nparts_,Bx,.true. )
    CALL GPIC_EulerToLag(this,this%lby_,this%nparts_,By,.false.)
    CALL GPIC_EulerToLag(this,this%lbz_,this%nparts_,Bz,.false.)
    CALL GPIC_EulerToLag(this,this%lfx_,this%nparts_,Ex,.false.)
    CALL GPIC_EulerToLag(this,this%lfy_,this%nparts_,Ey,.false.)
    CALL GPIC_EulerToLag(this,this%lfz_,this%nparts_,Ez,.false.)
!   Lorentz force
!$omp parallel do
    DO j = 1, this%nparts_
       this%lfx_(j) = this%lfx_(j)+this%pvy_(j)*this%lbz_(j)-this%pvz_(j)*this%lby_(j)
       this%lfy_(j) = this%lfy_(j)+this%pvz_(j)*this%lbx_(j)-this%pvx_(j)*this%lbz_(j)
       this%lfz_(j) = this%lfz_(j)+this%pvx_(j)*this%lby_(j)-this%pvy_(j)*this%lbx_(j)
    ENDDO

    ! ... x:
    dtfact = dt*xk*this%invdel_(1)
!$omp parallel do
    DO j = 1, this%nparts_
      this%px_ (j) = this%ptmp0_(1,j) + dtfact*this%pvx_(j)
      this%pvx_(j) = this%ttmp0_(1,j) + dtv*this%lfx_(j)
    ENDDO

    ! ... y:
    dtfact = dt*xk*this%invdel_(2)
!$omp parallel do
    DO j = 1, this%nparts_
      this%py_ (j) = this%ptmp0_(2,j) + dtfact*this%pvy_(j)
      this%pvy_(j) = this%ttmp0_(2,j) + dtv*this%lfy_(j)
    ENDDO

    ! ... z:
    dtfact = dt*xk*this%invdel_(3)
!$omp parallel do
    DO j = 1, this%nparts_
      this%pz_ (j) = this%ptmp0_(3,j) + dtfact*this%pvz_(j)
      this%pvz_(j) = this%ttmp0_(3,j) + dtv*this%lfz_(j)
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
!   DEALLOCATE(lid,gid)

    CALL ChargPIC_EndStageRKK(this,Ex,Ey,Ez,xk)

  END SUBROUTINE ChargPIC_StepRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE ChargPIC_EndStageRKK(this,vx,vy,vz,xk)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : EndStageRKK
!  DESCRIPTION: Called at the end of all RK-like stages to
!               complete inertial particle update.

!  ARGUMENTS  :
!    this    : 'this' class instance
!    Ez,Ey,Ez: compoments of velocity field, in real space, partially
!              updated, possibly. These will be overwritten!
!    xk      : multiplicative RK time stage factor
!-----------------------------------------------------------------
    USE fprecision
    USE commtypes
    USE mpivars
    USE grid

    IMPLICIT NONE
    CLASS(ChargPIC)    ,INTENT(INOUT)                      :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: vx,vy,vz
    REAL(KIND=GP),INTENT   (IN)                            :: xk
    INTEGER                                                :: j,ng

    ! u(t+dt) = u*: done already

    ! IF using nearest-neighbor interface, do particle exchange
    ! between nearest-neighbor tasks BEFORE z-PERIODIZING particle coordinates.
    ! Note this interface has not been tested yet for test particles.
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_NN ) THEN
      IF (this%nprocs_.GT.1) THEN
        CALL GTStart(this%htimers_(GPTIME_COMM))
        CALL this%gpcomm_%IdentifyExchV(this%id_,this%pz_,this%nparts_,ng,    &
                                       this%lxbnds_(3,1),this%lxbnds_(3,2)) 
        IF (ng.GT.this%partbuff_) THEN
          PRINT *, 'Rank', this%myrank_, 'resizing: nparts=', ng, ' | partbuff=',&
                   this%partbuff_, ' --> ', this%partbuff_ + &
                   (1+(ng-this%partbuff_)/this%partchunksize_)*this%partchunksize_
          this%partbuff_ = this%partbuff_ + &
                (1+(ng-this%partbuff_)/this%partchunksize_)*this%partchunksize_
          CALL this%ResizeArrays(this%partbuff_,.true.)
        END IF
        CALL MPI_BARRIER(this%comm_, this%ierr_)
        CALL this%gpcomm_%PartExchangeV(this%id_,this%px_,this%py_,this%pz_,  &
             this%nparts_,this%lxbnds_(3,1),this%lxbnds_(3,2),GPEXCH_INIT)
        CALL MPI_BARRIER(this%comm_, this%ierr_)
        CALL this%gpcomm_%PartExchangeV(this%id_,this%ptmp0_(1,:),            &
             this%ptmp0_(2,:),this%ptmp0_(3,:),this%nparts_,this%lxbnds_(3,1),&
             this%lxbnds_(3,2),GPEXCH_UPDT)
        CALL MPI_BARRIER(this%comm_, this%ierr_)
        CALL this%gpcomm_%PartExchangeV(this%id_,this%ttmp0_(1,:),            &
             this%ttmp0_(2,:),this%ttmp0_(3,:),this%nparts_,                  &
             this%lxbnds_(3,1),this%lxbnds_(3,2),GPEXCH_UPDT)
        CALL MPI_BARRIER(this%comm_, this%ierr_)
        CALL this%gpcomm_%PartExchangeV(this%id_,this%pvx_,this%pvy_,this%pvz_,&
             this%nparts_,this%lxbnds_(3,1),this%lxbnds_(3,2),GPEXCH_UPDT)
        CALL MPI_BARRIER(this%comm_, this%ierr_)
        CALL this%gpcomm_%PartExchangeV(this%id_,this%weight_,this%prop_, &
             this%prop_,this%nparts_,this%lxbnds_(3,1),this%lxbnds_(3,2), &
             GPEXCH_END)
        CALL GTAcc(this%htimers_(GPTIME_COMM))
      END IF
      ! Enforce periodicity in x and y:
      CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,3)
      ! Enforce periodicity in z and ptmp0(3):
      CALL GPart_MakePeriodicZ(this,this%pz_,this%ptmp0_(3,:),this%nparts_)
    ENDIF

    ! If using VDB interface, do synch-up, and get local work:
    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN

      ! Enforce periodicity in x, y, & z:
      CALL GPart_MakePeriodicP(this,this%px_,this%py_,this%pz_,this%nparts_,7)

      IF ( .NOT.GPart_PartNumConsistent(this,this%nparts_) ) THEN
        IF ( this%myrank_.eq.0 ) THEN
          WRITE(*,*) 'ChargPIC_EndStepRKK: Inconsistent particle count'
        ENDIF
      ENDIF
      ! Synch up VDB, if necessary:
      CALL GTStart(this%htimers_(GPTIME_COMM))
      CALL this%gpcomm_%VDBSynch(this%vdb_,this%maxparts_,this%id_, &
                     this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
      CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
                     this%pvx_,this%pvy_,this%pvz_,this%nparts_,this%ptmp1_)
      CALL GPart_CopyLocalWrk(this,this%pvx_,this%pvy_,this%pvz_, &
                     this%vdb_,this%gptmp0_,this%maxparts_)
      CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
                     this%ttmp0_(1,:),this%ttmp0_(2,:),this%ttmp0_(3,:),&
                     this%nparts_,this%ptmp1_)
      CALL GPart_CopyLocalWrk(this,this%ttmp0_(1,:),this%ttmp0_(2,:), &
                     this%ttmp0_(3,:),this%vdb_,this%gptmp0_,this%maxparts_)

      CALL this%gpcomm_%VDBSynch(this%ptmp0_,this%maxparts_,this%id_, &
               this%weight_,this%weight_,this%weight_,this%nparts_,this%ptmp1_)
      CALL GPIC_CopyLocalWrkScalar(this,this%weight_,this%vdb_,this%ptmp0_(1,:),this%maxparts_)
 
      CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
                     this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:),&
                     this%nparts_,this%ptmp1_)
      CALL GPart_GetLocalWrk_aux(this,this%id_,this%px_,this%py_,this%pz_,&
                       this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:),&
                       this%nparts_,this%vdb_,this%gptmp0_,this%maxparts_)
 
     CALL GTAcc(this%htimers_(GPTIME_COMM))

      CALL MPI_ALLREDUCE(this%nparts_,ng,1,MPI_INTEGER,   &
                         MPI_SUM,this%comm_,this%ierr_)

    ! vx, vy, vz were not used so far. Can be used in the future
    ! to compute acceleration, or as temporary arrays.

      IF ( this%myrank_.EQ.0 .AND. ng.NE.this%maxparts_) THEN
        WRITE(*,*)'ChargPIC_EndStepRKK: inconsistent d.b.: expected: ', &
                 this%maxparts_, '; found: ',ng
        STOP
      ENDIF

    ENDIF

    RETURN

  END SUBROUTINE ChargPIC_EndStageRKK
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE ChargPIC_InitFromFields(this,n,ux,uy,uz,T)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : InitFromFields
!  DESCRIPTION: Initializes particle position in uniform lattice 
!               and velocities with a Gaussian distribution 
!               function, with given mean local values and thermal 
!               speed given by local temperature. Other
!               parameters are initialized with GPart_Init.
!  ARGUMENTS:
!    this    : 'this' class instance
!    n       : particle density field to set each particle 
!              weight (not implemented yet)
!    ux,uy,uz: mean particle velocity field
!    T       : particle temperature field
!-----------------------------------------------------------------
    USE random
    USE var
    USE mpivars
    USE grid

    IMPLICIT NONE
    CLASS(ChargPIC)   , INTENT(INOUT)                   :: this
    REAL(KIND=GP),INTENT(IN),DIMENSION(nx,ny,ksta:kend) :: n,ux,uy,uz,T
    REAL(KIND=GP)                            :: gauss,vr,vth,w
    REAL(KIND=GP)                            :: del,delx,dely,delz
    DOUBLE PRECISION                         :: vmx,vmy,vmz
    INTEGER                                  :: ppc,pps,ib,lag,rppc
    INTEGER                                  :: i,j,k,ii,jj,kk
    INTEGER                                  :: ppx,ppy,ppz,d

    d = 3  ! Dimension of the system
    IF (nx.EQ.1) d = d - 1
    IF (ny.EQ.1) d = d - 1
    IF (nz.EQ.1) d = d - 1
    IF (d.EQ.0) THEN
      IF ( this%myrank_.eq.0 ) THEN
        WRITE(*,*) 'GPIC_InitFromFields: Number of cells must be greater than 1'
        STOP
      ENDIF
    END IF

    ppc = this%maxparts_/(nx*ny*nz)

    pps = ppc**(1.0/d)
    rppc = pps**d    ! Round to closest perfect d-power (trivial for d=1)

    this%nparts_ = nx*ny*(kend - ksta + 1)*ppc
    IF (this%nparts_.GT.this%partbuff_) THEN
      this%partbuff_ = this%partbuff_+(1+(this%nparts_-this%partbuff_) &
                              /this%partchunksize_)*this%partchunksize_
      CALL ChargPIC_ResizeArrays(this,this%partbuff_,.true.)
    END IF
    ib = nx*ny*(ksta-1)*ppc - 1
    lag = 1

    del = 1.0_GP/pps
    IF (nx.EQ.1) THEN
      ppx  = 1
      delx = 1.0_GP
    ELSE
      ppx  = pps
      delx = del
    END IF
    IF (ny.EQ.1) THEN
      ppy  = 1
      dely = 1.0_GP
    ELSE
      ppy  = pps
      dely = del
    END IF
    IF (nz.EQ.1) THEN
      ppz  = 1
      delz = 1.0_GP
    ELSE
      ppz  = pps
      delz = del
    END IF

    DO i = 1,nx
      DO j = 1,ny
        DO k = ksta,kend
          vmx = 0.0D0
          vmy = 0.0D0
          vmz = 0.0D0
          vth = SQRT(T(i,j,k))
          w   = n(i,j,k)*this%icv_
          DO ii = 1,ppx
            DO jj = 1,ppy
              DO kk = 1,ppz
                this%id_(lag) = lag + ib
                this%px_(lag) = (i-1.00_GP) + (ii-0.50_GP)*delx
                this%py_(lag) = (j-1.00_GP) + (jj-0.50_GP)*dely
                this%pz_(lag) = (k-1.00_GP) + (kk-0.50_GP)*delz
                this%weight_(lag) = w
                CALL random_gaussian(gauss)
                vr = gauss*vth
                this%pvx_(lag) = vr
                vmx = vmx + vr 
                CALL random_gaussian(gauss)
                vr = gauss*vth
                this%pvy_(lag) = vr
                vmy = vmy + vr 
                CALL random_gaussian(gauss)
                vr = gauss*vth
                this%pvz_(lag) = vr
                vmz = vmz + vr 
                lag = lag + 1
              END DO
            END DO
          END DO
          DO ii = rppc+1,ppc
            this%id_(lag) = lag + ib
            CALL prandom_number(this%px_(lag))
            CALL prandom_number(this%py_(lag))
            CALL prandom_number(this%pz_(lag))
            this%px_(lag) = this%px_(lag) + (i-1.00_GP)
            this%py_(lag) = this%py_(lag) + (j-1.00_GP)
            this%pz_(lag) = this%pz_(lag) + (k-1.00_GP)
            this%weight_(lag) = w
            CALL random_gaussian(gauss)
            vr = gauss*vth
            this%pvx_(lag) = vr
            vmx = vmx + vr
            CALL random_gaussian(gauss)
            vr = gauss*vth
            this%pvy_(lag) = vr
            vmy = vmy + vr
            CALL random_gaussian(gauss)
            vr = gauss*vth
            this%pvz_(lag) = vr
            vmz = vmz + vr
            lag = lag + 1
          END DO
          vmx = vmx/ppc - ux(i,j,k)
          vmy = vmy/ppc - uy(i,j,k)
          vmz = vmz/ppc - uz(i,j,k)
          lag = lag - ppc
          DO ii = 1,ppc
            this%pvx_(lag) = this%pvx_(lag) - vmx
            this%pvy_(lag) = this%pvy_(lag) - vmy
            this%pvz_(lag) = this%pvz_(lag) - vmz
            lag = lag + 1
          END DO
        END DO
      END DO
    END DO

    IF ( this%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN
       CALL this%gpcomm_%VDBSynch(this%vdb_,this%maxparts_,this%id_, &
                          this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
       CALL this%gpcomm_%VDBSynch(this%gptmp0_,this%maxparts_,this%id_, &
                          this%px_,this%py_,this%pz_,this%nparts_,this%ptmp1_)
       CALL GPart_GetLocalWrk(this,this%id_,this%px_,this%py_,this%pz_, &
                           this%nparts_,this%vdb_,this%maxparts_)

       CALL this%gpcomm_%VDBSynch(this%ptmp0_,this%maxparts_,this%id_, &
                          this%pvx_,this%pvy_,this%pvz_,this%nparts_,this%ptmp1_)
       CALL GPart_CopyLocalWrk(this,this%pvx_,this%pvy_,this%pvz_, &
                            this%vdb_,this%ptmp0_,this%maxparts_)

       CALL this%gpcomm_%VDBSynch(this%ptmp0_,this%maxparts_,this%id_, &
                          this%weight_,this%weight_,this%weight_,this%nparts_,this%ptmp1_)
       CALL GPIC_CopyLocalWrkScalar(this,this%weight_,this%vdb_,this%ptmp0_(1,:),this%maxparts_)

       IF ( this%wrtunit_ .EQ. 1 ) THEN ! rescale coordinates to box units
          this%ptmp0_(1,:) = this%vdb_(1,:)*this%delta_(1)
          this%ptmp0_(2,:) = this%vdb_(2,:)*this%delta_(2)
          this%ptmp0_(3,:) = this%vdb_(3,:)*this%delta_(3)
          CALL GPart_ascii_write_lag(this,1,'.','xlgInitRndSeed','000',0.0_GP,&
               this%maxparts_,this%ptmp0_(1,:),this%ptmp0_(2,:),this%ptmp0_(3,:))
       ELSE
          CALL GPart_ascii_write_lag(this,1,'.','xlgInitRndSeed','000',0.0_GP,&
               this%maxparts_,this%vdb_(1,:),this%vdb_(2,:),this%vdb_(3,:))
       ENDIF

       IF ( .NOT.GPart_PartNumConsistent(this,this%nparts_) ) THEN
         IF ( this%myrank_.eq.0 ) THEN
           WRITE(*,*) 'GPIC_InitLattice: Invalid particle after GetLocalWrk call'
           STOP
         ENDIF
       ENDIF
    END IF

    RETURN

  END SUBROUTINE ChargPIC_InitFromFields
!-----------------------------------------------------------------
!-----------------------------------------------------------------
  
  SUBROUTINE random_gaussian(gauss)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : random_gaussian
!  DESCRIPTION: Generates random number following gaussian 
!               distribution.
!  ARGUMENTS:
!    gauss  : 
!-----------------------------------------------------------------
    USE random
    USE var

    IMPLICIT NONE
    REAL(KIND=GP),INTENT(OUT)        :: gauss
    REAL(KIND=GP)                    :: u1,u2,twopi,low

    twopi = 2.0_GP*pi
    low   = 1.0e-20_GP

    CALL prandom_number(u1)
    CALL prandom_number(u2)
    u1 =  MAX( u1, low)
    u2 =  twopi * u2
    gauss =  sqrt( -2.0_GP * log( u1)) * CMPLX( cos(u2), sin(u2))

  END SUBROUTINE random_gaussian
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE ChargPIC_GetTemperature(this,T)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetTemperature
!  DESCRIPTION: Deposits particle kinetic energy density into grid.
!
!  ARGUMENTS  :
!    this     : 'this' class instance
!    T        : Eulerian field containing particle temperature (OUT)
!-----------------------------------------------------------------
    USE grid
    USE mpivars

    IMPLICIT NONE
    CLASS(ChargPIC),INTENT(INOUT)                          :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: T
    INTEGER                                                :: lag

    DO lag=1,this%nparts_
      this%prop_(lag) = this%weight_(lag)*(this%pvx_(lag)**2&
                        +this%pvy_(lag)**2+this%pvz_(lag)**2)
    END DO

    CALL GPIC_LagToEuler(this,this%prop_,this%nparts_,T,.true.)

    RETURN

  END SUBROUTINE ChargPIC_GetTemperature
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE ChargPIC_GetTemperatureAnis(this,Tpa,Tpe,dir)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : GetTemperature
!  DESCRIPTION: Deposits particle parallel and perpendicular 
!               kinetic energy density into grid.
!
!  ARGUMENTS  :
!    this     : 'this' class instance
!    Tpa,Tpe  : Eulerian field containing particle parallel and
!               perpendicular temperature (OUT)
!    dir      : Component defining 'parallel' (0->x, 1->y, 2->z)
!-----------------------------------------------------------------
    USE grid
    USE mpivars

    IMPLICIT NONE
    CLASS(ChargPIC),INTENT(INOUT)                          :: this
    REAL(KIND=GP),INTENT(INOUT),DIMENSION(nx,ny,ksta:kend) :: Tpa,Tpe
    INTEGER      ,INTENT(IN)                               :: dir
    INTEGER                                                :: lag

    IF (dir.EQ.0) THEN
      DO lag=1,this%nparts_
        this%prop_(lag) = this%weight_(lag)*(this%pvy_(lag)**2&
                                            +this%pvz_(lag)**2)
      END DO
      CALL GPIC_LagToEuler(this,this%prop_,this%nparts_,Tpe,.false.)

      DO lag=1,this%nparts_
        this%prop_(lag) = this%weight_(lag)*this%pvx_(lag)**2
      END DO
      CALL GPIC_LagToEuler(this,this%prop_,this%nparts_,Tpa,.false.)
    ELSE IF (dir.EQ.1) THEN
      DO lag=1,this%nparts_
        this%prop_(lag) = this%weight_(lag)*(this%pvx_(lag)**2&
                                            +this%pvz_(lag)**2)
      END DO
      CALL GPIC_LagToEuler(this,this%prop_,this%nparts_,Tpe,.false.)

      DO lag=1,this%nparts_
        this%prop_(lag) = this%weight_(lag)*this%pvy_(lag)**2
      END DO
      CALL GPIC_LagToEuler(this,this%prop_,this%nparts_,Tpa,.false.)
    ELSE IF (dir.EQ.2) THEN
       DO lag=1,this%nparts_
        this%prop_(lag) = this%weight_(lag)*(this%pvx_(lag)**2&
                                            +this%pvy_(lag)**2)
      END DO
      CALL GPIC_LagToEuler(this,this%prop_,this%nparts_,Tpe,.false.)

      DO lag=1,this%nparts_
        this%prop_(lag) = this%weight_(lag)*this%pvz_(lag)**2
      END DO
      CALL GPIC_LagToEuler(this,this%prop_,this%nparts_,Tpa,.false.)
    END IF

    RETURN

  END SUBROUTINE ChargPIC_GetTemperatureAnis
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  SUBROUTINE ChargPIC_ResizeArrays(this,new_size,onlyinc,exc)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  METHOD     : Resize_Arrays
!  DESCRIPTION: Resize all arrays in the ChargPIC class (including
!               subclases, i.e. communicator, spline)
!  ARGUMENTS  :
!    this    : 'this' class instance
!    new_size: new number of particles
!    onlyinc : if true, will only resize to increase array size
!-----------------------------------------------------------------
!$  USE threads

    IMPLICIT NONE
    CLASS(ChargPIC) ,INTENT(INOUT)                      :: this
    INTEGER         ,INTENT(IN)                         :: new_size
    LOGICAL         ,INTENT(IN)                         :: onlyinc
    LOGICAL         ,INTENT(IN)   ,OPTIONAL             :: exc
    INTEGER                                             :: n

    CALL VGPIC_ResizeArrays(this,new_size,onlyinc,exc)

    n = SIZE(this%lfx_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank1(this%lfx_,new_size,.false.)
    END IF
    n = SIZE(this%lfy_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank1(this%lfy_,new_size,.false.)
    END IF
    n = SIZE(this%lfz_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank1(this%lfz_,new_size,.false.)
    END IF

    n = SIZE(this%lbx_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank1(this%lbx_,new_size,.false.)
    END IF
    n = SIZE(this%lby_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank1(this%lby_,new_size,.false.)
    END IF
    n = SIZE(this%lbz_)
    IF ((n.lt.new_size).OR.((n.gt.new_size).AND..NOT.onlyinc)) THEN
      CALL Resize_ArrayRank1(this%lbz_,new_size,.false.)
    END IF

    RETURN

  END SUBROUTINE ChargPIC_ResizeArrays
!-----------------------------------------------------------------
!-----------------------------------------------------------------
